#ifndef SGLBM_H
#define SGLBM_H
#include <../../src/uq2D.h>
#include <../../src/uq2D.hh>
#include <iostream>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <sys/stat.h>
#include <iomanip>

#define totalCollocation
// #define onlyUVCollocation
#define stochastic_omega
// #define constant_omega

template<typename T>
class sglbm{
public:
  // ---------------- physical & lattice parameters ------------------------
  int  N{0}, nx{0}, ny{0};
  T    L{1}, dx{1}, lx{0}, ly{0};
  T    physVelocity{0}, physViscosity{0}, physDensity{1};
  T    Re{0}, tau{0}, dt{0};
  T    conversionViscosity{1}, conversionVelocity{1}, conversionDensity{1},
        conversionMass{1}, conversionForce{1}, u0{0}, omega0{0};
  std::string dir;

  // ---------------- stochastic data --------------------------------------
  std::unique_ptr<GeneralizedPolynomialChaos<T>> ops;   // shared pattern
  int  No{0}, total_nq{0};
  std::vector<int> polynomial_types;
  std::vector<T>   parameter1, parameter2;

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
    const Parameters&   params,
    UncertaintyQuantification<T>& uq)
    : dir(_dir) {

    ops        = uq.getOps();
    No         = ops->getPolynomialsOrder();
    total_nq   = ops->getQuadraturePointsNumber();
    parameter1 = params.parameter1;
    parameter2 = params.parameter2;
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

  std::vector<T> find_intersection(std::vector<T> center, T radius, std::vector<int> start_point, std::vector<int> end_point)
  {
    T x0 = center[0];
    T y0 = center[1];

    T x1 = start_point[0];
    T y1 = start_point[1];

    T x2 = end_point[0];
    T y2 = end_point[1];

    std::vector<T> inp(2, 0.0);
    if (radius == 0) {
      return { x1, y1 };
    }

    if (x1 == x2) {
      inp[0] = x1;
      if (std::abs(radius) >= std::abs(x1 - x0)) {
        T p1 = y0 - std::sqrt(radius * radius - (x1 - x0) * (x1 - x0));
        T p2 = y0 + std::sqrt(radius * radius - (x1 - x0) * (x1 - x0));
        if (std::max(y1, y2) >= p2) {
          inp[2] = p2;
        }
        if (std::min(y1, y2) <= p1) {
          inp[2] = p1;
        }
      }
    }
    else {
      T k = (y1 - y2) / (x1 - x2);
      T b0 = y1 - k * x1;

      T a = k * k + 1;
      T b = 2.0 * k * (b0 - y0) - 2.0 * x0;
      T c = (b0 - y0) * (b0 - y0) + x0 * x0 - radius * radius;
      T delta = b * b - 4 * a * c;
      if (delta >= 0) {
        T p1x = (-1.0 * b - std::sqrt(delta)) / (2 * a);
        T p2x = (-1.0 * b + std::sqrt(delta)) / (2 * a);
        T p1y = k * p1x + b0;
        T p2y = k * p2x + b0;
        if (p1x >= std::min(x1, x2) && p1x <= std::max(x1, x2)) {
          inp[1] = p1x;
          inp[2] = p1y;
        }
        else {
          inp[1] = p2x;
          inp[2] = p2y;
        }
      }
    }
    return inp;
  }


  void UnitConverterFromResolutionAndRelaxationTime(const Parameters params)
  {
    std::cout << "set simulation parameter" << std::endl;
    N = params.resolution;
    tau = params.tau;
    L = params.L;
    physVelocity = params.physVelocity;

    if (params.Re == 0) {
      physViscosity = params.physViscosity;
      Re = physVelocity * L / physViscosity;
    }
    else if (params.physViscosity == 0){
      Re = params.Re;
      physViscosity = physVelocity * L / Re;
    }

    physDensity = 1.0;
    
    dx = L / N;
    dt = (tau - 0.5) / 3.0 * (dx * dx) / physViscosity;

    conversionViscosity = dx * dx / dt;
    conversionVelocity = dx / dt;
    conversionDensity = 1.0;
    conversionMass = conversionDensity * dx * dx * dx;
    conversionForce = conversionMass * dx / dt / dt;
    u0 = physVelocity / conversionVelocity;
    omega0 = 1.0 / (3 * physViscosity / conversionViscosity + 0.5);

    std::cout << "resolution " << N << std::endl;
    std::cout << "tau = " << tau << std::endl;
    std::cout << "nu = " << physViscosity << std::endl;
    std::cout << "Re = " << Re << std::endl;
    std::cout << "u0 = " << u0 << std::endl;
    std::cout << "converstionVelocity = " << conversionVelocity << std::endl;
    std::cout << "converstionViscosity = " << conversionViscosity << std::endl;
    std::cout << "dx = " << dx << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "Ma = " << u0 * std::sqrt(3.0) << std::endl;
  }

  void setCircle(T centerX, T centerY, T radius)
  {
    bouzidiQ.resize(nx);
    for (int i = 0; i < nx; ++i) {
      bouzidiQ[i].resize(ny);
      for (int j = 0; j < ny; ++j) {
        bouzidiQ[i][j].resize(9);
      }
    }

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        if (material[i][j] == 3) {
          for (int k = 1; k < 9; ++k) {
            std::vector<T> intersection(2, 0.0);
            intersection = find_intersection({ centerX / dx + 1,centerY / dx + 1 }, L / dx, { i,j }, { i + cx[k],j + cy[k] });
            bouzidiQ[i][j][k] = 1.0 - std::sqrt((i - intersection[1]) * (i - intersection[1]) + (j - intersection[2]) * (j - intersection[2])) / sqrt(cx[k] * cx[k] + cy[k] * cy[k]);
          }
        }
      }
    }
  }

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

    std::vector<T> rRan(total_nq, physDensity),
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

        outputFile << x << "\t" << y << "\t" << std::sqrt(ops->mean(mag)) << "\t" << ops->mean(u[id]) * conversionVelocity << "\t" << ops->mean(v[id]) * conversionVelocity << "\t" << std::sqrt(ops->std(mag)) << "\t" << ops->std(u[id]) * conversionVelocity << "\t" << ops->std(v[id]) * conversionVelocity  << "\t" << material[i][j] <<  "\n";
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

  T getParameter1(int index) {
    return parameter1[index];
  }

  T getParameter2(int index) {
    return parameter2[index];
  }

  int get_polynomials_order() {
    return No;
  }

  int get_quadrature_points_number() {
    return total_nq;
  }
  
};



#endif // LBM_H
