#ifndef SGLBM_H
#define SGLBM_H
#include "polynomials.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <sys/stat.h>
#include <iomanip>

#define totalCollocation
// #define onlyUVCollocation

class sglbm{
public:
  int N;
  double L;
  double dx;
  double lx;
  double ly;
  int nx;
  int ny;
  double physVelocity;
  double physViscosity;
  double physDensity;
  double Re;
  double tau;
  double dt;
  double conversionViscosity;
  double conversionVelocity;
  double conversionDensity;
  double conversionMass;
  double conversionForce;
  double u0;
  double omega0;
  std::string dir;

  polynomials ops;


  double cs2 = 1.0 / 3.0;
  std::vector<int> cx = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
  std::vector<int> cy = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
  std::vector<double> w = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
  std::vector<int> kinv = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

  std::vector<double> conversionViscosityChaos;
  std::vector<double> conversionVelocityChaos;
  std::vector<double> conversionDensityChaos;
  std::vector<double> conversionMassChaos;
  std::vector<double> conversionForceChaos;
  std::vector<double> ReChaos;
  std::vector<double> u0Chaos;
  std::vector<double> dtChaos;
  std::vector<double> tauChaos;
  std::vector<double> omegaChaos;

  std::vector<double> conversionViscosityRan;
  std::vector<double> conversionVelocityRan;
  std::vector<double> conversionDensityRan;
  std::vector<double> conversionMassRan;
  std::vector<double> conversionForceRan;
  std::vector<double> ReRan;
  std::vector<double> u0Ran;
  std::vector<double> dtRan;
  std::vector<double> tauRan;

  

  std::vector<std::vector<int>> material;
  std::vector<std::vector<std::vector<double>>> bouzidiQ;
  std::vector<std::vector<std::vector<double>>> u;
  std::vector<std::vector<std::vector<double>>> v;
  std::vector<std::vector<std::vector<double>>> rho;
  std::vector<std::vector<std::vector<std::vector<double>>>> f;
  std::vector<std::vector<std::vector<std::vector<double>>>> F;
  std::vector<std::vector<std::vector<std::vector<double>>>> feq;


  sglbm(std::string _dir, const Parameters& params)
        : ops(params.nq, params.order, params.parameter1, params.parameter2, params.polynomialType, params.points_weights_method) 
  {
    dir = _dir;

  }

  std::vector<double> find_intersection(std::vector<double> center, double radius, std::vector<int> start_point, std::vector<int> end_point)
  {
    double x0 = center[0];
    double y0 = center[1];

    double x1 = start_point[0];
    double y1 = start_point[1];

    double x2 = end_point[0];
    double y2 = end_point[1];

    std::vector<double> inp(2, 0.0);
    if (radius == 0) {
      return { x1, y1 };
    }

    if (x1 == x2) {
      inp[0] = x1;
      if (std::abs(radius) >= std::abs(x1 - x0)) {
        double p1 = y0 - std::sqrt(radius * radius - (x1 - x0) * (x1 - x0));
        double p2 = y0 + std::sqrt(radius * radius - (x1 - x0) * (x1 - x0));
        if (std::max(y1, y2) >= p2) {
          inp[2] = p2;
        }
        if (std::min(y1, y2) <= p1) {
          inp[2] = p1;
        }
      }
    }
    else {
      double k = (y1 - y2) / (x1 - x2);
      double b0 = y1 - k * x1;

      double a = k * k + 1;
      double b = 2.0 * k * (b0 - y0) - 2.0 * x0;
      double c = (b0 - y0) * (b0 - y0) + x0 * x0 - radius * radius;
      double delta = b * b - 4 * a * c;
      if (delta >= 0) {
        double p1x = (-1.0 * b - std::sqrt(delta)) / (2 * a);
        double p2x = (-1.0 * b + std::sqrt(delta)) / (2 * a);
        double p1y = k * p1x + b0;
        double p2y = k * p2x + b0;
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


  void UnitConverterFromResolutionAndRelaxationTime(Parameters params)
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

  void setCircle(double centerX, double centerY, double radius)
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
            std::vector<double> intersection(2, 0.0);
            intersection = find_intersection({ centerX / dx + 1,centerY / dx + 1 }, L / dx, { i,j }, { i + cx[k],j + cy[k] });
            bouzidiQ[i][j][k] = 1.0 - std::sqrt((i - intersection[1]) * (i - intersection[1]) + (j - intersection[2]) * (j - intersection[2])) / sqrt(cx[k] * cx[k] + cy[k] * cy[k]);
          }
        }
      }
    }
  }

  void prepareLattice()
  {
    std::cout << "starting initialization" << std::endl;

    u.resize(nx);
    v.resize(nx);
    rho.resize(nx);
    f.resize(nx);
    F.resize(nx);
    feq.resize(nx);

    for (int i = 0; i < nx; ++i) {
      u[i].resize(ny);
      v[i].resize(ny);
      rho[i].resize(ny);
      f[i].resize(ny);
      F[i].resize(ny);
      feq[i].resize(ny);
      for (int j = 0; j < ny; ++j) {
        u[i][j].resize(ops.No);
        v[i][j].resize(ops.No);
        rho[i][j].resize(ops.No);
        
        for (int alpha = 0; alpha < ops.No; ++alpha) {
          u[i][j][alpha] = 0.0;
          v[i][j][alpha] = 0.0;
          rho[i][j][alpha] = 0.0;
        }

        f[i][j].resize(9);
        F[i][j].resize(9);
        feq[i][j].resize(9);
        for (int k = 0; k < 9; ++k) {
          f[i][j][k].resize(ops.No);
          F[i][j][k].resize(ops.No);
          feq[i][j][k].resize(ops.No);
          for (int alpha = 0; alpha < ops.No; ++alpha) {
            f[i][j][k][alpha] = 0.0;
            F[i][j][k][alpha] = 0.0;
            feq[i][j][k][alpha] = 0.0;
          }
        }
      }
    }
    
    omegaChaos.resize(ops.No);
  }

  void initializeDistributionFunction() {
    
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          feq[i][j][k][0] = equilibrium(rho[i][j][0], u[i][j][0], v[i][j][0], cx[k], cy[k], w[k]);
          f[i][j][k][0] = feq[i][j][k][0];
          F[i][j][k][0] = feq[i][j][k][0];
        }
      }
    }
    
  }

  void collision()
  {

    std::vector<double> QSlice(ops.No, 0.0);
    //std::cout << "collision start" << std::endl;
#pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          collisionTerm(f[i][j][k], feq[i][j][k], omegaChaos, QSlice);

          for (int alpha = 0; alpha < ops.No; ++alpha) {
            F[i][j][k][alpha] = f[i][j][k][alpha] + QSlice[alpha];
          }
        }
      }
    }
  }

  void collisionTerm(std::vector<double> _f, std::vector<double> _feq, std::vector<double> _omega, std::vector<double>& Q)
  {

    for (int i = 0; i < ops.No; ++i) {
      double sum1 = 0.0;
      double sum2 = 0.0;

      for (int j = 0; j < ops.No; ++j) {
        for (int k = 0; k < ops.No; ++k) {
          size_t flatIndex = i + ops.No * (k + ops.No * j);
          sum1 += _omega[j] * _feq[k] * ops.t3Product[flatIndex];
          sum2 += _omega[j] * _f[k] * ops.t3Product[flatIndex];
        }
      }

      Q[i] = (sum1 - sum2) * ops.t2Product_inv[i];
    }
  }

  double equilibrium(double _r, double _u, double _v, int _cx, int _cy, double _w)
  {
    double t1 = _u * _u + _v * _v;
    double t2 = _u * _cx + _v * _cy;
    return _r * _w * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
  }

  void streaming()
  {
    //std::vector<std::vector<std::vector<std::vector<double>>>> ftmp(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(9, std::vector<double>(order + 1, 0.0))));
    //std::cout << "streaming start" << std::endl;
    #pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          int ii = (i + nx + cx[k]) % (nx);
          int jj = (j + ny + cy[k]) % (ny);

          f[ii][jj][k] = F[i][j][k];
        }
      }
    }

  }


  void reconstruction()
  {
    //std::cout<<"max threads: " << nProcessors<<std::endl;
    //omp_set_dynamic(0);     // Explicitly disable dynamic teams
    std::vector<double> rRan(ops.total_nq, 0.0);
    std::vector<double> uRan(ops.total_nq, 0.0);
    std::vector<double> vRan(ops.total_nq, 0.0);
    std::vector<double> ruRan(ops.total_nq, 0.0);
    std::vector<double> rvRan(ops.total_nq, 0.0);

    std::vector<double> fSlice(ops.No, 0.0);
    std::vector<double> feqSlice(ops.No, 0.0);

#pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        std::vector<double> rChaos(ops.No, 0.0);
        std::vector<double> ruChaos(ops.No, 0.0);
        std::vector<double> rvChaos(ops.No, 0.0);
        std::vector<double> uChaos(ops.No, 0.0);
        std::vector<double> vChaos(ops.No, 0.0);

        std::vector<double> rRan(ops.total_nq, 0.0);
        std::vector<double> uRan(ops.total_nq, 0.0);
        std::vector<double> vRan(ops.total_nq, 0.0);
        std::vector<double> ruRan(ops.total_nq, 0.0);
        std::vector<double> rvRan(ops.total_nq, 0.0);

        for (int alpha = 0; alpha < ops.No; ++alpha) {
          for (int k = 0; k < 9; ++k) {
            rChaos[alpha] += f[i][j][k][alpha];
          }
        }
        rho[i][j] = rChaos;

        ops.chaos2ran(rChaos, rRan);

        if (material[i][j] == 1) {
          for (int alpha = 0; alpha < ops.No; ++alpha) {
            for (int k = 0; k < 9; ++k) {
              ruChaos[alpha] += f[i][j][k][alpha] * cx[k];
              rvChaos[alpha] += f[i][j][k][alpha] * cy[k];
            }
          }            
          //std::cout << "check" << std::endl;
          ops.chaos2ran(ruChaos, ruRan);
          ops.chaos2ran(rvChaos, rvRan);

          for (int sample = 0; sample < ops.total_nq; sample++) {
            uRan[sample] = ruRan[sample] / rRan[sample];
            vRan[sample] = rvRan[sample] / rRan[sample];
          }

          ops.ran2chaos(uRan, uChaos);
          ops.ran2chaos(vRan, vChaos);

          u[i][j] = uChaos;
          v[i][j] = vChaos;
        }
        else if (material[i][j] == 2)
        {
          for (int alpha = 0; alpha < ops.No; ++alpha) {
            u[i][j][alpha] = 0.0;
            v[i][j][alpha] = 0.0;
          }     
          
          for (int sample = 0; sample < ops.total_nq; ++sample) {
            uRan[sample] = 0.0;
            vRan[sample] = 0.0;
          }
        }
        // Stochastic velocity boundary condition
        else if (material[i][j] == 3){
          ops.chaos2ran(u[i][j], uRan);
          std::vector<double> empty(ops.total_nq, 0.0);
          vRan.swap(empty);
        }

        for (int k = 0; k < 9; ++k) {
          #if defined(onlyUVCollocation)
            std::vector<double> ruuChaos(ops.No, 0.0);
            std::vector<double> rvvChaos(ops.No, 0.0);
            std::vector<double> ruvChaos(ops.No, 0.0);

            ops.chaos_product(ruChaos, uChaos, ruuChaos);
            ops.chaos_product(rvChaos, vChaos, rvvChaos);
            ops.chaos_product(ruChaos, vChaos, ruvChaos);

            for (int alpha = 0; alpha < ops.No; ++alpha) {
              double t1 = ruuChaos[alpha] + rvvChaos[alpha];
              double t2 = ruChaos[alpha] * cx[k] + rvChaos[alpha] * cy[k];
              double t22 = ruuChaos[alpha] * cx[k] * cx[k] + 2 * ruvChaos[alpha] * cx[k] * cy[k] + rvvChaos[alpha] * cy[k] * cy[k];
              feq[i][j][k][alpha] = w[k] * (rChaos[alpha] + 3 * t2 + 4.5 * t22 - 1.5 * t1);
            }
          #elif defined(totalCollocation)
            std::vector<double> feqRan(ops.total_nq, 0.0);
            for (int sample = 0; sample < ops.total_nq; sample++) {
              feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
            }

            ops.ran2chaos(feqRan, feqSlice);
            // for (int alpha = 0; alpha < ops.No; ++alpha) {
            //   feq[i][j][k][alpha] = feqSlice[alpha];
            // }
            feq[i][j][k] = feqSlice;
            feqRan.clear();
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

    outputFile << "variables = \"X\", \"Y\", \"RhoMean\", \"uMean\", \"vMean\", \"RhoStd\", \"uStd\", \"vStd\", \"geometry\"\n";
    outputFile << "ZONE I = " << nx << ", J = " << ny << ", F = POINT\n";

    std::vector<double> rSlice(ops.No, 0.0);
    std::vector<double> uSlice(ops.No, 0.0);
    std::vector<double> vSlice(ops.No, 0.0);

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x = i * dx;
        double y = j * dx;

        for (int alpha = 0; alpha < ops.No; ++alpha) {
          rSlice[alpha] = rho[i][j][alpha];
          uSlice[alpha] = u[i][j][alpha];
          vSlice[alpha] = v[i][j][alpha];
        }
        outputFile << x << "\t" << y << "\t" << ops.mean(rho[i][j]) << "\t" << ops.mean(u[i][j]) * conversionVelocity << "\t" << ops.mean(v[i][j]) * conversionVelocity << "\t" << ops.std(rSlice) << "\t" << ops.std(uSlice) * conversionVelocity << "\t" << ops.std(vSlice) * conversionVelocity  << "\t" << material[i][j] <<  "\n";
      }
    }
    outputFile.close();
  }
    

    void boundary()
    {
      std::vector<double> rRan(ops.total_nq, 0.0);
      std::vector<double> uRan(ops.total_nq, 0.0);
      std::vector<double> vRan(ops.total_nq, 0.0);
      std::vector<double> rSlice(ops.No, 0.0);
      std::vector<double> uSlice(ops.No, 0.0);
      std::vector<double> feqSlice(ops.No, 0.0);
      std::vector<double> chaos(2,0.0);
      ops.convert2affinePCE(ops.parameter1[0]*u0, ops.parameter2[0]*u0, ops.polynomial_types[0],chaos);
      uSlice[0] = chaos[0];
      uSlice[1] = chaos[1];
      ops.chaos2ran(uSlice, uRan);

      #pragma omp for collapse(2)
      for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
          if (material[i][j] == 2){
            for (int k = 0; k < 9; ++k){
              int new_i = i+cx[k];
              int new_j = j+cy[k];
              if ((new_i!=-1) && (new_i!=nx) && (new_j!=0) && (new_j!=ny)){
                if (material[new_i][new_j] == 1){
                  for (int alpha = 0; alpha < ops.No; ++alpha){
                    F[i][j][k][alpha] = F[new_i][new_j][kinv[k]][alpha];
                  }
                }
              }
            }
          }
          else if (material[i][j] == 3){
            for (int alpha = 0; alpha < ops.No; ++alpha) {
              rSlice[alpha] = rho[i][j-1][alpha];
            }
            ops.chaos2ran(rSlice, rRan);

            for (int k = 0; k < 9; ++k) {
              std::vector<double> feqRan(ops.total_nq, 0.0);
              for (int sample = 0; sample < ops.total_nq; sample++) {
                feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
              }

              ops.ran2chaos(feqRan, feqSlice);
              for (int alpha = 0; alpha < ops.No; ++alpha) {
                F[i][j][k][alpha] = feqSlice[alpha] + F[i][j-1][k][alpha] - feq[i][j-1][k][alpha];
              }
            }
          }
        }
      }
    }

};



#endif // LBM_H
