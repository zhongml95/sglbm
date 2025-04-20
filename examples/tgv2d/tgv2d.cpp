// -----------------------------------------------------------------------------
// Taylor–Green vortex 2‑D driver (updated for the refactored sglbm)
// -----------------------------------------------------------------------------
//   * Works with the flattened SoA layout inside sglbm
//   * Uses the new convenience accessors   u_at(i,j) / v_at(i,j) / rho_at(i,j)
//   * No per‑node heap allocations in the main loop
//   * Cleaned‑up "auto dist" shadowing and minor compilation warnings
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <filesystem>


#include "../../src/sglbm.h"
#include "../../src/postprocessing_sglbm.h"

using T = double;

// -----------------------------------------------------------------------------
//  Compile‑time switches for the stochastic parameter you want to study
// -----------------------------------------------------------------------------
// #define stochastic_Re
#define stochastic_viscosity
// #define stochastic_velocity

// -----------------------------------------------------------------------------
//  Helpers
// -----------------------------------------------------------------------------

template<typename T>
T calc_tke_error(const sglbm<T>& sglbm, int count)
{
    std::vector<T> tke(sglbm.ops->getPolynomialsOrder(), 0.0);
    T tkeAna = 0.0;
    totalKineticEnergy(sglbm, tke, tkeAna, count);
    return std::abs((sglbm.ops->mean(tke) - tkeAna) / tkeAna);
}

void setGeometry(sglbm<T>& s)
{
    std::cout << "start setting geometry\n";

    s.nx = s.N;
    s.ny = s.N;

    /* create an nx × ny array filled with 1 (fluid) */
    s.material.assign(s.nx, std::vector<int>(s.ny, 1));

    std::cout << "nx: " << s.nx << ", ny: " << s.ny << '\n';
    std::cout << "finish setting geometry\n";
}


// -----------------------------------------------------------------------------
//  Initial condition
// -----------------------------------------------------------------------------

template<typename T>
void initialize(sglbm<T>& s, [[maybe_unused]] const Parameters& params)
{
    s.prepareLattice();

    std::cout << "polynomial order: " << s.get_polynomials_order() << "\n";
    std::cout << "quadrature points number: " << s.get_quadrature_points_number() << "\n";

    //--------------------------------------------------------------------------
    //  Set relaxation parameter omega (possibly stochastic)
    //--------------------------------------------------------------------------
    {
        std::vector<T> omegaRan(s.get_quadrature_points_number());
        std::vector<T> omegaChaos(s.get_polynomials_order(), 0.0);

    #if defined(stochastic_Re)
        // Re ~ U[Re*par1, Re*par2]
        std::vector<T> ReChaos(s.get_polynomials_order(), 0.0);
        std::vector<T> ReRan(s.get_quadrature_points_number());
        std::array<T,2> coeff{};
        auto distRe = uniform(s.Re * params.parameter1[0], s.Re * params.parameter2[0]);
        s.ops->convert2affinePCE(distRe, coeff);
        ReChaos[0] = coeff[0];
        ReChaos[1] = coeff[1];
        s.ops->chaosToRandom(ReChaos, ReRan);
        for (int q = 0; q < s.get_quadrature_points_number(); ++q)
            omegaRan[q] = 1.0 / (3.0 * (s.physVelocity * s.L / ReRan[q]) / s.conversionViscosity + 0.5);

    #elif defined(stochastic_viscosity)
        // viscosity provided as SC samples from the GPC object
        auto physViscRan = s.ops->getStochasticCollocationSample(); // returns [nq][?]
        for (int q = 0; q < s.get_quadrature_points_number(); ++q)
            omegaRan[q] = 1.0 / (3.0 * physViscRan[q][0] / s.conversionViscosity + 0.5);

    #else   // deterministic omega (default)
        omegaRan.assign(s.get_quadrature_points_number(),
                        1.0 / (3.0 * (s.physVelocity * s.L / s.Re) / s.conversionViscosity + 0.5));
    #endif

        s.ops->randomToChaos(omegaRan, omegaChaos);
        s.omegaChaos = omegaChaos;

        std::cout << "mean omega: " << s.ops->mean(omegaChaos)
                  << ", std(omega): " << s.ops->std(omegaChaos) << "\n";
    }

    //--------------------------------------------------------------------------
    //  Taylor–Green vortex initial fields
    //--------------------------------------------------------------------------
    for (int i = 0; i < s.nx; ++i) {
        for (int j = 0; j < s.ny; ++j) {
            const T x = i * s.dx;
            const T y = j * s.dx;

            std::vector<T> uChaos(s.get_polynomials_order(), 0.0);
            std::vector<T> vChaos(s.get_polynomials_order(), 0.0);
            std::vector<T> rChaos(s.get_polynomials_order(), 0.0);

            // deterministic part (0‑th chaos coefficient)
            rChaos[0] = 1.0 - 1.5 * s.u0 * s.u0 * std::cos(x + y) * std::cos(x - y);
            uChaos[0] = -s.u0 * std::cos(x) * std::sin(y);
            vChaos[0] =  s.u0 * std::sin(x) * std::cos(y);

            // store in lattice
            s.rho_at(i,j) = rChaos;
            s.u_at(i,j)   = uChaos;
            s.v_at(i,j)   = vChaos;

        #if defined(stochastic_velocity)
            // add velocity perturbations on chaos modes 1‑4
            auto addPerturbation = [&](int mode, T a, T b) {
                std::array<T,2> coeff{};
                auto dist = uniform(a,b);
                s.ops->convert2affinePCE(dist, coeff);
                s.u_at(i,j)[mode] -= s.u0 * 0.25 * coeff[1] * std::cos(x) * std::sin(y);
                s.v_at(i,j)[mode] += s.u0 * 0.25 * coeff[1] * std::sin(x) * std::cos(y);
            };

            addPerturbation(1, params.parameter1[0], params.parameter2[0]);
            addPerturbation(2, params.parameter1[1], params.parameter2[1]);
            addPerturbation(3, params.parameter1[2], params.parameter2[2]);
            addPerturbation(4, params.parameter1[3], params.parameter2[3]);
        #endif
        }
    }

    s.initializeDistributionFunction();
    std::cout << "finish initializing" << std::endl;
}

// -----------------------------------------------------------------------------
//  Main
// -----------------------------------------------------------------------------

int main([[maybe_unused]] int  argc,
  [[maybe_unused]] char* argv[])
{
    Parameters params;
    readParameters("./parameters.dat", params);

    // Recompute tau from user‑supplied viscosity / Re
    const T physViscosity = params.physVelocity * params.L / params.Re;
    params.tau = 3 * physViscosity + 0.5;

    // data directory ---------------------------------------------------------
    const std::string dir     = "./data/tgv/t5/ReNr" + std::to_string(params.order) +
                                "Nq" + std::to_string(params.nq) +
                                "N"  + std::to_string(params.resolution) + "/";
    const std::string dirAna  = dir + "final/";

    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    std::filesystem::create_directories(dirAna);

    //--------------------------------------------------------------------------
    //  Uncertainty‑quantification setup (GPC)
    //--------------------------------------------------------------------------
    UncertaintyQuantification<T> uq(UQMethod::GPC);
    auto distVisc = uniform(params.parameter1[0] * physViscosity,
                            params.parameter2[0] * physViscosity);
    uq.initializeGPC(params.order, params.nq, distVisc);

    //--------------------------------------------------------------------------
    //  Build LBM object
    //--------------------------------------------------------------------------
    sglbm<T> s(dir, params, uq);
    s.UnitConverterFromResolutionAndRelaxationTime(params);

    setGeometry(s);
    initialize(s, params);

    //--------------------------------------------------------------------------
    //  Time stepping
    //--------------------------------------------------------------------------
    const T td = 1.0 / (s.physViscosity * (s.dx * s.dx * 2.0));
    const int maxIter = static_cast<int>(td * 0.5);
    const int outputInterval = 100;

    size_t cores = omp_get_num_procs();
    omp_set_dynamic(0);
    omp_set_num_threads(static_cast<int>(cores));
    std::cout << "num Threads: " << cores << std::endl;

    s.output(dir, 0);

    double tStart = omp_get_wtime();
    double tCheckpoint = tStart;
    int    iterOut = 0;

#pragma omp parallel
    for (int iter = 1; iter <= maxIter; ++iter) {
        s.collision();
        s.streaming();
        s.reconstruction();

#pragma omp single
        {
            if (iter % outputInterval == 0 || iter == maxIter) {
                double now = omp_get_wtime();
                T err = calc_tke_error(s, iter);
                std::cout << "iter: " << iter
                          << ", Δt: " << now - tCheckpoint << " s"
                          << ", TKE err: " << err << std::endl;
                s.output(dir, iter);
                tCheckpoint = now;
                iterOut = iter;
            }
        }
    }

    // final output -----------------------------------------------------------                
    double now = omp_get_wtime();
                T err = calc_tke_error(s, iterOut);
                std::cout << "iter: " << iterOut
                          << ", Δt: " << now - tCheckpoint << " s"
                          << ", TKE err: " << err << std::endl;
    s.output(dir, iterOut);
    double tEnd = omp_get_wtime();
    T errFinal = calc_tke_error(s, iterOut);

    std::cout << "total time: " << tEnd - tStart << " s, final TKE err: " << errFinal << std::endl;

    velocity_central(dir, s, (s.nx/2)+1, (s.ny/2)+1);
    velocity_all(dir, s);
    outputTKE(dir, s, iterOut, tEnd - tStart);

    return 0;
}
