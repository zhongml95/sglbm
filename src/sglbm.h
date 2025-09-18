#pragma once
#ifndef SGLBM_H
#define SGLBM_H

#include "uq.h"
#include "uq.hh"

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "parameters.h"
#include "unitconverter.h"



#define totalCollocation
// #define onlyUVCollocation
#define stochastic_omega
// #define constant_omega

using namespace olb;
using namespace olb::uq;

template<typename T>
class sglbm{
public:
  // ---------------- physical & lattice parameters ------------------------
  int  N{0}, nx{0}, ny{0};
  T    L{1}, dx{1}, lx{0}, ly{0};
  std::string dir;

  UnitConverter<T> unit;

  // ---------------- stochastic data --------------------------------------
  std::unique_ptr<GeneralizedPolynomialChaos<T>> ops;   // shared pattern
  int  No{0}, total_nq{0};
  std::vector<int> polynomial_types;

  // ---------------- discrete velocity model ------------------------------
  const T cs2{ T(1) / 3 };                                   // c_s^2
  const std::array<int, 9> cx{ 0, 1, 0, -1, 0, 1, -1, -1,  1 };
  const std::array<int, 9> cy{ 0, 0, 1,  0,-1, 1,  1, -1, -1 };
  const std::array<int, 9> kinv{ 0, 3, 4, 1, 2, 7, 8, 5, 6 };
  const std::array<T,   9> w{ 4./9, 1./9, 1./9, 1./9, 1./9,
                              1./36,1./36,1./36,1./36 };

  // ---------------- chaos scalars (stochastic) ---------------------------
  std::vector<T> omegaChaos;

  // ---------------- geometry and fields ----------------------------------
  std::vector<std::vector<int>> material;            // nx × ny (Eulerian mask)
  std::vector<std::vector<std::vector<T>>> bouzidiQ; // curved walls (rare)

  /* flattened (cells = nx*ny) ------------------------------------------- */
  std::vector<std::vector<T>> u, v, rho;                    // cells × No
  std::vector<std::vector<std::vector<T>>> f, F, feq;       // cells × 9 × No


  // ======================================================================
  //  Ctor
  // ======================================================================
  sglbm(const std::string& _dir,
        const UnitConverter<T>& unitConverter,
        UncertaintyQuantification<T>& uq)
  : dir(_dir), unit(unitConverter)
  {
    ops = std::make_unique<olb::uq::GeneralizedPolynomialChaos<T>>(uq.getOps());
    No         = ops->getPolynomialsOrder();
    total_nq   = ops->getQuadraturePointsNumber();

    N = unit.getResolution();
    nx = unit.getResolution(); // update it in the future if needed
    ny = unit.getResolution();
  }

  // ======================================================================
  //  Helpers
  // ======================================================================
  inline size_t idx(int i,int j) const { return static_cast<size_t>(i)*ny + j; }

  /* compatibility wrappers ------------------------------------------------*/
  inline std::vector<T>&       u_at  (int i,int j){ return  u[idx(i,j)]; }
  inline const std::vector<T>& u_at  (int i,int j) const{ return  u[idx(i,j)]; }
  inline std::vector<T>&       v_at  (int i,int j){ return  v[idx(i,j)]; }
  inline const std::vector<T>& v_at  (int i,int j) const{ return  v[idx(i,j)]; }
  inline std::vector<T>&       rho_at(int i,int j){ return rho[idx(i,j)]; }
  inline const std::vector<T>& rho_at(int i,int j) const{ return rho[idx(i,j)]; }

    // ----------------------------------------------------------------------
    // Field allocation (SoA)
    // ----------------------------------------------------------------------
    void prepareLattice()
    {
        const size_t cells = static_cast<size_t>(nx)*ny;
        u.assign  (cells, std::vector<T>(No, 0));
        v.assign  (cells, std::vector<T>(No, 0));
        rho.assign(cells, std::vector<T>(No, 0));

        f .resize(cells);
        F .resize(cells);
        feq.resize(cells);
        for (size_t id = 0; id < cells; ++id) {
            f  [id].assign(9, std::vector<T>(No, 0));
            F  [id] = f[id];
            feq[id] = f[id];
        }
        omegaChaos.assign(No, 0.0);
    }

    // ----------------------------------------------------------------------
    // Initial equilibrium distribution   (thread‑safe)
    // ----------------------------------------------------------------------
  void initializeDistributionFunction() {

    std::vector<T> rRan(total_nq, unit.getPhysDensity()),
    uRan(total_nq, 0.0), vRan(total_nq, 0.0),
    feqRan(total_nq), feqChaos(No);
    
    for (size_t id = 0; id < u.size(); ++id) {
      ops->chaosToRandom(u[id], uRan);
      ops->chaosToRandom(v[id], vRan);

      for (int k = 0; k < 9; ++k) {
          for (int q = 0; q < total_nq; ++q)
              feqRan[q] = equilibrium(rRan[q], uRan[q], vRan[q],
                                       cx[k], cy[k], w[k]);
          ops->randomToChaos(feqRan, feqChaos);
          feq[id][k] = f[id][k] = F[id][k] = feqChaos;
      }
    }
  }

      // ----------------------------------------------------------------------
    //  Collision (BGK; stochastic omega)  — thread‑local GPC
    // ----------------------------------------------------------------------
  void collision()
  {
    std::vector<T> Q(No, 0.0);
    auto& gpc = *ops;            // shared but read‑only for chaos ops
    #pragma omp for schedule(static)
    for (size_t id = 0; id < f.size(); ++id) {;
        for (int k = 0; k < 9; ++k) {
          gpc.chaosProduct(omegaChaos,MatrixOperations::vectorSubtraction(feq[id][k], f[id][k]),Q);
          F[id][k] = MatrixOperations::vectorAdd(f[id][k], Q);
        }
      }
  }

    // ----------------------------------------------------------------------
    //  Streaming (periodic)
    // ----------------------------------------------------------------------  
  void streaming()
  {
    // std::cout << "streaming start" << std::endl;
    #pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        const size_t idSrc = idx(i,j);
        for (int k = 0; k < 9; ++k) {
          int ii = (i + nx + cx[k]) % (nx);
          int jj = (j + ny + cy[k]) % (ny);

          f[idx(ii,jj)][k] = F[idSrc][k];
        }
      }
    }

  }

  T equilibrium(T _r, T _u, T _v, int _cx, int _cy, T _w)
  {
    T t1 = _u * _u + _v * _v;
    T t2 = _u * _cx + _v * _cy;
    return _r * _w * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
  }

  std::vector<T> equilibrium(std::vector<T> _r, std::vector<T> _u, std::vector<T> _v, int _cx, int _cy, T _w)
  {
    std::vector<T> _feq(No, 0.0);
    std::vector<T> _u2(No, 0.0);
    std::vector<T> _v2(No, 0.0);
    std::vector<T> _feq_without_rho(No, 0.0);
    std::vector<T> _t2(No, 0.0);
    ops->chaosProduct(_u, _u, _u2);
    ops->chaosProduct(_v, _v, _v2);

    for (int i = 0; i < No; ++i) {
      _t2[i] = _u[i] * _cx + _v[i] * _cy;
    }
    std::vector<T> _t2_2(No, 0.0);
    ops->chaosProduct(_t2, _t2, _t2_2);

    for (int i = 0; i < No; ++i) {
      T t1 = _u2[i] + _v2[i];
      _feq_without_rho[i] = _w * (1.0 + 3.0 * _t2[i] + 4.5 * _t2_2[i] - 1.5 * t1);
    }

    ops->chaosProduct(_r, _feq_without_rho, _feq);
    return _feq;
  }


  bool check_positivity() {
    
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int alpha = 0; alpha < No; ++alpha) {
          if (rho[i][j][alpha] < 0) {
            std::cout << "alpha: " << alpha << ", rho: "<<  rho[i][j][alpha] << std::endl;
            // return false;
          }
        }
      }
    }
    return true;
  }

  void reconstruction()
  {
    /* thread‑local scratch buffers ------------------------------ */
    std::vector<T> rChaos(No), ruChaos(No), rvChaos(No);
    std::vector<T> uChaos(No), vChaos(No), feqSlice(No);

    std::vector<T> rRan(total_nq), uRan(total_nq),
                    vRan(total_nq), ruRan(total_nq), rvRan(total_nq),
                    feqRan(total_nq);


  #pragma omp for collapse(2)
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      const size_t id = idx(i,j);
      std::fill(rChaos.begin(), rChaos.end(), T(0));
      for (int alpha = 0; alpha < No; ++alpha) {
        for (int k = 0; k < 9; ++k) {
          rChaos[alpha] += f[id][k][alpha];
        }
      }
      rho[id] = rChaos;
      ops->chaosToRandom(rChaos, rRan);

      if (material[i][j] == 1) {                  // fluid
        std::fill(ruChaos.begin(), ruChaos.end(), T(0));
        std::fill(rvChaos.begin(), rvChaos.end(), T(0));

        for (int alpha = 0; alpha < No; ++alpha) {
          for (int k = 0; k < 9; ++k) {
            ruChaos[alpha] += f[id][k][alpha] * cx[k];
            rvChaos[alpha] += f[id][k][alpha] * cy[k];
          }
        }
        ops->chaosToRandom(ruChaos, ruRan);
        ops->chaosToRandom(rvChaos, rvRan);

        for (int sample = 0; sample < total_nq; sample++) {
          uRan[sample] = ruRan[sample] / rRan[sample];
          vRan[sample] = rvRan[sample] / rRan[sample];
        }

        ops->randomToChaos(uRan, uChaos);
        ops->randomToChaos(vRan, vChaos);

        u[id] = uChaos;
        v[id] = vChaos;
      }
      else if (material[i][j] == 2) {             // solid
        std::fill(u[id].begin(), u[id].end(), T(0));
        std::fill(v[id].begin(), v[id].end(), T(0));
        std::fill(uRan.begin(), uRan.end(), T(0));
        std::fill(vRan.begin(), vRan.end(), T(0));
      }
      else if (material[i][j] == 3) {             // inlet
        ops->chaosToRandom(u[id], uRan);
        std::fill(vRan.begin(), vRan.end(), T(0));
      }

  /* -------- build local equilibrium ------------------- */
      for (int k = 0; k < 9; ++k) {
        #if defined(onlyUVCollocation)
          std::vector<T> ruuChaos(No, 0.0);
          std::vector<T> rvvChaos(No, 0.0);
          std::vector<T> ruvChaos(No, 0.0);

          ops->chaosProduct(ruChaos, uChaos, ruuChaos);
          ops->chaosProduct(rvChaos, vChaos, rvvChaos);
          ops->chaosProduct(ruChaos, vChaos, ruvChaos);

          for (int alpha = 0; alpha < No; ++alpha) {
            T t1 = ruuChaos[alpha] + rvvChaos[alpha];
            T t2 = ruChaos[alpha] * cx[k] + rvChaos[alpha] * cy[k];
            T t22 = ruuChaos[alpha] * cx[k] * cx[k] + 2 * ruvChaos[alpha] * cx[k] * cy[k] + rvvChaos[alpha] * cy[k] * cy[k];
            feq[id][k][alpha] = w[k] * (rChaos[alpha] + 3 * t2 + 4.5 * t22 - 1.5 * t1);
          }
        #elif defined(totalCollocation)
          for (int sample = 0; sample < total_nq; sample++) {
            feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
          }

          ops->randomToChaos(feqRan, feqSlice);
          feq[id][k] = feqSlice;
        #endif
        }
      }
    }
  }

  void output(std::string dir, int iter)
  {
    std::cout << "output iter " << iter << " start" << std::endl;
    std::string filename = dir + std::to_string(iter) + ".dat";
    std::ofstream outputFile(filename);
    if (!outputFile) {
      std::cerr << "Error opening the file: " << filename << std::endl;
      return;
    }

    outputFile << "variables = \"X\", \"Y\", \"magMean\", \"uMean\", \"vMean\", \"magStd\", \"uStd\", \"vStd\", \"geometry\"\n";
    outputFile << "ZONE I = " << nx << ", J = " << ny << ", F = POINT\n";

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        const size_t id = idx(i,j);
        T x = i * dx;
        T y = j * dx;
        std::vector<T> u2(No, 0.0);
        std::vector<T> v2(No, 0.0);
        std::vector<T> mag(No, 0.0);
        ops->chaosProduct(u[id], u[id], u2);
        ops->chaosProduct(v[id], v[id], v2);
        ops->chaosSum(u2, v2, mag);

        outputFile << x << "\t" << y << "\t" << std::sqrt(ops->mean(mag)) << "\t" << ops->mean(u[id]) * unit.getConversionVelocity() << "\t" << ops->mean(v[id]) * unit.getConversionVelocity() << "\t" << std::sqrt(ops->std(mag)) << "\t" << ops->std(u[id]) * unit.getConversionVelocity() << "\t" << ops->std(v[id]) * unit.getConversionVelocity() << "\t" << material[i][j] <<  "\n";
      }
    }
    outputFile.close();
  }

  void boundary()
  {
    std::vector<T> rRan(total_nq, 0.0);
    std::vector<T> uRan(total_nq, 0.0);
    std::vector<T> vRan(total_nq, 0.0);
    std::vector<T> rSlice(No, 0.0);
    std::vector<T> uSlice(No, 0.0);
    std::vector<T> feqSlice(No, 0.0);
    std::vector<T> chaos(2,0.0);
    auto dist = ops->getDistributions();
    ops->convert2affinePCE(dist[0], chaos);
    uSlice[0] = chaos[0];
    uSlice[1] = chaos[1];
    ops->chaosToRandom(uSlice, uRan);

    #pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i){
      for (int j = 0; j < ny; ++j){
        const size_t id = idx(i,j);
        if (material[i][j] == 2){
          for (int k = 0; k < 9; ++k){
            int new_i = i+cx[k];
            int new_j = j+cy[k];
            if ((new_i!=-1) && (new_i!=nx) && (new_j!=0) && (new_j!=ny)){
              if (material[new_i][new_j] == 1){
                for (int alpha = 0; alpha < No; ++alpha){
                  F[idx][k][alpha] = F[idx(new_i, new_j)][kinv[k]][alpha];
                }
              }
            }
          }
        }
        else if (material[i][j] == 3){
          for (int alpha = 0; alpha < No; ++alpha) {
            rSlice[alpha] = rho[i][j-1][alpha];
          }
          ops->chaosToRandom(rSlice, rRan);

          for (int k = 0; k < 9; ++k) {
            std::vector<T> feqRan(total_nq, 0.0);
            for (int sample = 0; sample < total_nq; sample++) {
              feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
            }

            ops->randomToChaos(feqRan, feqSlice);
            for (int alpha = 0; alpha < No; ++alpha) {
              F[id][k][alpha] = feqSlice[alpha] + F[idx(i, j-1)][k][alpha] - feq[idx(i, j-1)][k][alpha];
            }
          }
        }
      }
    }
  }

  int get_polynomials_order() {
    return No;
  }

  int get_quadrature_points_number() {
    return total_nq;
  }

};



#endif // SGLBM_H
