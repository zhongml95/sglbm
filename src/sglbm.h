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
  double dy;
  double lx;
  double ly;
  int nx;
  int ny;
  double physVelocity;
  double Re;
  double physViscosity;
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
  std::string exm;

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
  std::vector<std::vector<std::vector<double>>> omega;
  std::vector<std::vector<std::vector<std::vector<double>>>> f;
  std::vector<std::vector<std::vector<std::vector<double>>>> F;
  std::vector<std::vector<std::vector<std::vector<double>>>> feq;


    sglbm(std::string _dir, std::string _exm, const Parameters& params)
        : ops(params.nq, params.order, params.parameter1, params.parameter2, params.polynomialType, params.points_weights_method) 
  {
    dir = _dir;
    exm = _exm;

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

  void setGeometry(double _L, double _N, double _lx, double _ly, const std::vector<std::vector<int>>& _material) {
    std::cout << "start setting geometry" << std::endl;
    L = _L;
    N = _N;
    dx = L / N;
    dy = L / N;
    lx = _lx;
    ly = _ly;
    
    nx = int(lx / dx) + 1;
    ny = int(ly / dy) + 1;

    if (exm == "tgv"){
      nx = int(lx / dx);
      ny = int(ly / dy);
    }

    material.resize(nx);
    for (int i = 0; i < nx; ++i) {
      material[i].resize(ny);
      for (int j = 0; j < ny; ++j) {
        material[i][j] = _material[i][j];
      }
    }
    std::cout << "finish setting geometry" << std::endl;

  }

  void setFluid(double _physVelocity, double _nu, double _tau)
  {
    std::cout << "set simulation parameter" << std::endl;
    physVelocity = _physVelocity;
    //Re = _Re;
    physViscosity = _nu;//physVelocity*L/Re;
    Re = physVelocity * L / physViscosity;
    tau = _tau;
    dt = (tau - 0.5) / 3.0 * (dx * dx) / physViscosity;
    conversionViscosity = dx * dx / dt;
    conversionVelocity = dx / dt;
    conversionDensity = 1.0;
    conversionMass = conversionDensity * dx * dx * dx;
    conversionForce = conversionMass * dx / dt / dt;
    u0 = physVelocity / conversionVelocity;
    omega0 = 1.0 / (3 * physViscosity / conversionViscosity + 0.5);
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

  void initialize()
  {
    std::cout << "starting initialization" << std::endl;

    std::cout << "resolution " << ny << std::endl;
    std::cout << "nx = " << nx << "\t" << "ny = " << ny << std::endl;
    std::cout << "tau = " << tau << std::endl;
    std::cout << "nu = " << physViscosity << std::endl;
    std::cout << "Re = " << Re << std::endl;
    std::cout << "u0 = " << u0 << std::endl;
    std::cout << "converstionVelocity = " << conversionVelocity << std::endl;
    std::cout << "converstionViscosity = " << conversionViscosity << std::endl;

    u.resize(nx);
    v.resize(nx);
    rho.resize(nx);
    omega.resize(nx);
    f.resize(nx);
    F.resize(nx);
    feq.resize(nx);

    for (int i = 0; i < nx; ++i) {
      u[i].resize(ny);
      v[i].resize(ny);
      rho[i].resize(ny);
      omega[i].resize(ny);
      f[i].resize(ny);
      F[i].resize(ny);
      feq[i].resize(ny);
      for (int j = 0; j < ny; ++j) {
        u[i][j].resize(ops.No);
        v[i][j].resize(ops.No);
        rho[i][j].resize(ops.No);
        omega[i][j].resize(ops.No);
        for (int alpha = 0; alpha < ops.No; ++alpha) {
          u[i][j][alpha] = 0.0;
          v[i][j][alpha] = 0.0;
          rho[i][j][alpha] = 0.0;
          omega[i][j][alpha] = 0.0;
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

    ReChaos.resize(ops.No);
    ReRan.resize(ops.total_nq);


    // std::vector<double> uChaos(ops.No, 0.0);
    // std::vector<double> vChaos(ops.No, 0.0);
    // std::vector<double> rChaos(ops.No, 0.0);
    std::vector<double> omegaChaos(ops.No, 0.0);
    std::vector<double> physViscosityChaos(ops.No, 0.0);
    std::vector<double> chaos(2, 0.0);
    //rChaos[0] = 1.0;
    std::cout << "omega: " << 1.0 / (3 * physViscosity / conversionViscosity + 0.5) << std::endl;

    if (exm == "tgv"){
      std::vector<double> omegaRan(ops.total_nq, 0.0);
      std::vector<double> physViscosityRan(ops.total_nq, 0.0);

      //if(ops.polynomialType == 0) {
      ops.convert2affinePCE(Re * ops.parameter1[0], Re * ops.parameter2[0], ops.polynomial_types[0], chaos);
      //}
      //else if (ops.polynomialType == 1) {
      //  ops.convert2affinePCE(physViscosity * ops.parameter1, physViscosity * ops.parameter2, chaos);
      //}

      ReChaos[0] = chaos[0];
      ReChaos[1] = chaos[1]; 

      std::cout << "ReChaos: ";
      for (int alpha = 0; alpha < ops.No; ++alpha) {
        std::cout << ReChaos[alpha] << "\t";
      }
      std::cout << std::endl;

      ops.chaos2ran(ReChaos, ReRan);


      for (int sample = 0; sample < ops.total_nq; ++sample) {
        physViscosityRan[sample] = physVelocity * L / ReRan[sample];
        omegaRan[sample] = 1.0 / ( 3.0 * physViscosityRan[sample] / conversionViscosity + 0.5 );
      }

      ops.ran2chaos(omegaRan, omegaChaos);

      //ops.convert2affinePCE(ops.parameter1 / ( 3.0 * physViscosity / conversionViscosity + 0.5 ), ops.parameter2 / ( 3.0 * physViscosity / conversionViscosity + 0.5 ), chaos);
      //omegaChaos[0] = chaos[0];
      //omegaChaos[1] = chaos[1];
      std::cout << "omegaChaos: ";
      for (int alpha = 0; alpha < ops.No; ++alpha) {
        std::cout << omegaChaos[alpha] << "\t";
      }
      std::cout << std::endl;
      /*if(ops.polynomialType == 0) {
        ops.convert2affinePCE(1.0 / (3 * physViscosity * ops.parameter2 / conversionViscosity + 0.5), 1.0 / (3 * physViscosity * ops.parameter1 / conversionViscosity + 0.5),chaos);
      }
      else if (ops.polynomialType == 1) {
        ops.convert2affinePCE((1.0 / (3 * physViscosity / conversionViscosity + 0.5)) * ops.parameter1, (1.0 / (3 * physViscosity / conversionViscosity + 0.5)) * ops.parameter2, chaos);
      }*/
    
      std::cout << "chaos: " << omegaChaos[0] << "\t" << omegaChaos[1] << "\t" << omegaChaos[2] << std::endl;
      // chaos.clear();
    }
    else if (exm == "cavity2d" || exm == "stgv") {
      omegaChaos[0] = 1.0 / ( 3.0 * physViscosity / conversionViscosity + 0.5 );
    }


    //omp_set_num_threads(4);
//#pragma omp parallel for collapse(2)
    //#pragma omp single
    //{
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        std::vector<double> uChaos(ops.No+1, 0.0);
        std::vector<double> vChaos(ops.No+1, 0.0);
        std::vector<double> rChaos(ops.No+1, 0.0);
        // std::cout << "i: " << i << "\t" << "j: " << j << ": " << vChaos[2] << std::endl;
        //TGV            
        for (int alpha = 0; alpha < ops.No; ++alpha) {
          //printf("i: %d, j: %d, alpha: %d, omp_get_thread_num: %d\n", i,j,alpha,omp_get_thread_num()); 
          if(exm == "tgv"){
              double x = i * dx;
              double y = j * dy;
              
              //rChaos[alpha] = 1.0 - 1.5 * u0Chaos[alpha] * u0Chaos[alpha] * std::cos(x + y) * std::cos(x - y);
              //uChaos[alpha] = -u0Chaos[alpha] * std::cos(x) * std::sin(y);
              //vChaos[alpha] = u0Chaos[alpha] * std::sin(x) * std::cos(y);

              rChaos[0] = 1.0 - 1.5 * u0 * u0 * std::cos(x + y) * std::cos(x - y);
              uChaos[0] = -u0 * std::cos(x) * std::sin(y);
              vChaos[0] = u0 * std::sin(x) * std::cos(y);

              //omegaChaos[0] = omega0;
              //printf("i: %d, j: %d, alpha: %d, omp_get_thread_num: %d\n", i,j,alpha,omp_get_thread_num()); 
          }
          else if (exm == "cavity2d"){
            rChaos[0] = 1.0;
            omegaChaos[0] = 1.0 / (3 * physViscosity / conversionViscosity + 0.5);
          }
          else if (exm == "stgv"){

              double x = i * dx;
              double y = j * dy;

              rChaos[0] = 1.0 - 1.5 * u0 * u0 * std::cos(x + y) * std::cos(x - y);
              uChaos[0] = -u0 * std::cos(x) * std::sin(y);
              vChaos[0] = u0 * std::sin(x) * std::cos(y);
              std::vector<double> u0Chaos(ops.random_number_dimension + 1, 0.0);

              if (alpha <= ops.random_number_dimension && alpha != 0) {
                int random_id = alpha - 1;
                // if(alpha == 3)
                // std::cout << ops.parameter1[2] << "\t" << ops.parameter2[2] << "\t" << ops.polynomial_types[2] << std::endl;
                ops.convert2affinePCE(ops.parameter1[random_id], ops.parameter2[random_id], ops.polynomial_types[random_id], chaos);
                if (random_id == 0) {
                  u0Chaos[0] += u0 * chaos[0] * std::sin(x) * std::sin(y);
                  u0Chaos[alpha] = u0 * chaos[1] * std::sin(x) * std::sin(y);
                }
                else if (random_id == 1)  {
                  u0Chaos[0] += u0 * chaos[0] * std::sin(x) * std::cos(y);
                  u0Chaos[alpha] = u0 * chaos[1] * std::sin(x) * std::cos(y);
                }
                else if (random_id == 2) {
                  u0Chaos[0] += u0 * chaos[0] * std::cos(x) * std::sin(y);
                  // std::cout << std::cos(x) * std::sin(y) << std::endl;
                  // std::cout << u0 * chaos[1] << std::endl;
                  // std::cout << u0 * chaos[1] * std::cos(x) * std::sin(y) << std::endl;
                  u0Chaos[alpha] = u0 * chaos[1] * std::cos(x) * std::sin(y);
                }
                else if (random_id == 3) {
                  u0Chaos[0] += u0 * chaos[0] * std::cos(x) * std::cos(y);
                  u0Chaos[alpha] = u0 * chaos[1] * std::cos(x) * std::cos(y);
                }
              }
              if (alpha < ops.random_number_dimension+1 && alpha != 0) {

                uChaos[alpha] += 0.25 * u0Chaos[alpha] * std::cos(x) * std::sin(y);
                vChaos[alpha] += 0.25 * u0Chaos[alpha] * std::sin(x) * std::cos(y);
              }

          }
          
          u[i][j][alpha] = uChaos[alpha];
          v[i][j][alpha] = vChaos[alpha];
          rho[i][j][alpha] = rChaos[alpha];
          omega[i][j][alpha] = omegaChaos[alpha];

          // if(alpha < ops.random_number_dimension+1) {
          //   std::cout << uChaos[alpha] << std::endl;
          //   std::cout << u[i][j][alpha] << "\t";
          // }
        }
        // std::cout << std::endl;
      }
    }


    std::vector<double> rRan(ops.total_nq, 0.0);
    std::vector<double> uRan(ops.total_nq, 0.0);
    std::vector<double> vRan(ops.total_nq, 0.0);

    //std::vector<double> rSlice(ops.No + 1, 0.0);
    //std::vector<double> uSlice(ops.No + 1, 0.0);
    //std::vector<double> vSlice(ops.No + 1, 0.0);

    std::vector<double> feqRan(ops.total_nq, 0.0);
    std::vector<double> feqSlice(ops.No, 0.0);

    //#pragma omp parallel for num_threads(4)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        /*for (int alpha = 0; alpha < ops.No + 1; ++alpha) {
          rSlice[alpha] = rho[i][j][alpha];
          uSlice[alpha] = u[i][j][alpha];
          vSlice[alpha] = v[i][j][alpha];
        }*/

        ops.chaos2ran(rho[i][j], rRan);
        ops.chaos2ran(u[i][j], uRan);
        ops.chaos2ran(v[i][j], vRan);

        for (int k = 0; k < 9; ++k) {

          for (int sample = 0; sample < ops.total_nq; sample++) {
            feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
          }

          ops.ran2chaos(feqRan, feqSlice);
          for (int alpha = 0; alpha < ops.No; ++alpha) {
            feq[i][j][k][alpha] = feqSlice[alpha];
            //std::cout << feqSlice[alpha] << std::endl;
            f[i][j][k][alpha] = feqSlice[alpha];
            F[i][j][k][alpha] = feqSlice[alpha];
          }
        }
      }
    }

    //std::cout << "initialize finished" << std::endl;

    // uChaos.clear();
    // vChaos.clear();
    // rChaos.clear();
    // omegaChaos.clear();

    // rRan.clear();
    // uRan.clear();
    // vRan.clear();
    // feqRan.clear();
    // rSlice.clear();
    // uSlice.clear();
    // vSlice.clear();
    // feqSlice.clear();
    // std::cout << "clear data" << std::endl;

  }

    void collision()
  {

    //std::vector<double> omegaSlice(ops.No + 1, 0.0);
    //std::vector<double> fSlice(ops.No + 1, 0.0);
    //std::vector<double> feqSlice(ops.No + 1, 0.0);
    std::vector<double> QSlice(ops.No, 0.0);
    //std::cout << "collision start" << std::endl;
#pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          /*for (int alpha = 0; alpha < ops.No + 1; ++alpha) {
            omegaSlice[alpha] = omega[i][j][alpha];
            fSlice[alpha] = f[i][j][k][alpha];
            feqSlice[alpha] = feq[i][j][k][alpha];
            //std::cout << feq[i][j][k][alpha] << std::endl;
          }*/
          collisionTerm(f[i][j][k], feq[i][j][k], omega[i][j], QSlice);

          for (int alpha = 0; alpha < ops.No; ++alpha) {
            F[i][j][k][alpha] = f[i][j][k][alpha] + QSlice[alpha];
            //std::cout << F[i][j][k][alpha] << std::endl;
          }
        }
      }
    }
    //std::cout << "collision finished" << std::endl;

    //omegaSlice.clear();
    //fSlice.clear();
    //feqSlice.clear();
    //QSlice.clear();
  }

  void collisionTerm(std::vector<double> _f, std::vector<double> _feq, std::vector<double> _omega, std::vector<double>& Q)
  {

    for (int i = 0; i < ops.No; ++i) {
      double sum1 = 0.0;
      double sum2 = 0.0;

      for (int j = 0; j < ops.No; ++j) {
        for (int k = 0; k < ops.No; ++k) {
          // sum1 += _omega[j] * _feq[k] * ops.t3Product[j][k][i];
          // sum2 += _omega[j] * _f[k] * ops.t3Product[j][k][i];
          size_t flatIndex = i + ops.No * (k + ops.No * j);
          sum1 += _omega[j] * _feq[k] * ops.t3Product[flatIndex];
          sum2 += _omega[j] * _f[k] * ops.t3Product[flatIndex];
        }
      }

      Q[i] = (sum1 - sum2) / ops.t2Product[i];
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
          for (int alpha = 0; alpha < ops.No; ++alpha) {
            f[ii][jj][k][alpha] = F[i][j][k][alpha];
          }
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
          rho[i][j][alpha] = rChaos[alpha];
        }
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

          for (int alpha = 0; alpha < ops.No; ++alpha) {
            u[i][j][alpha] = uChaos[alpha];
            v[i][j][alpha] = vChaos[alpha];
          }
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
        else if (material[i][j] == 3){
          std::vector<double> chaos(2,0.0);
          ops.convert2affinePCE(ops.parameter1[0]*u0, ops.parameter2[0]*u0, ops.polynomial_types[0],chaos);
          uChaos[0] = chaos[0];
          uChaos[1] = chaos[1];
          ops.chaos2ran(uChaos, uRan);
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
            for (int alpha = 0; alpha < ops.No; ++alpha) {
              feq[i][j][k][alpha] = feqSlice[alpha];
            }
            feqRan.clear();
          #endif
        }
      }     
    }

  }

 void output(std::string dir, int iter, double total_computational_time)
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
        double y = j * dy;

        for (int alpha = 0; alpha < ops.No; ++alpha) {
          rSlice[alpha] = rho[i][j][alpha];
          uSlice[alpha] = u[i][j][alpha];
          vSlice[alpha] = v[i][j][alpha];
        }
        outputFile << x << "\t" << y << "\t" << ops.mean(rho[i][j]) << "\t" << ops.mean(u[i][j]) * conversionVelocity << "\t" << ops.mean(v[i][j]) * conversionVelocity << "\t" << ops.std(rSlice) << "\t" << ops.std(uSlice) * conversionVelocity << "\t" << ops.std(vSlice) * conversionVelocity  << "\t" << material[i][j] <<  "\n";

        // rSlice.clear();
        // uSlice.clear();
        // vSlice.clear();
      }
    }
    outputFile.close();

    if (exm == "tgv" || exm == "stgv"){
      std::string filenameTKE = dir + "final/tke.dat";
      std::ofstream outputFileTKE(filenameTKE);
      std::vector<double> tke(ops.No,0.0);
      double tkeAna = 0.0;
      totalKineticEnergy(tke, tkeAna, iter+1);

      outputFileTKE.precision(20);  
      outputFileTKE << ops.mean(tke) << "\t" << ops.std(tke) << "\t" << tkeAna << "\t" << total_computational_time;
      outputFileTKE.close();
      // tke.clear();
    }

    std::string filenameU = dir + "final/u.dat";
    std::ofstream outputFileU(filenameU);

    for(int j = 0; j < ny; ++j){
      int idx = (nx+1)/2;
      if (exm == "tgv" || exm == "stgv"){
        idx = nx / 2 + 1;
      }
      //std::vector<double> uSlice(ops.No, 0.0);
      for (int alpha = 0; alpha < ops.No; ++alpha) {
        uSlice[alpha] = u[idx][j][alpha];
      }

      outputFileU.precision(20);  
      outputFileU << j * dy << "\t" << ops.mean(uSlice) * conversionVelocity << "\t" << ops.std(uSlice) * conversionVelocity;

      if (exm == "tgv" || exm == "stgv"){
        double x = idx * dx;
        double y = j * dy;
        double k2 = dx*dx + dy*dy;
        double damp = std::exp(-k2*physViscosity*iter);
        outputFileU << "\t" << -u0 * std::cos(x) * std::sin(y) * damp * conversionVelocity;
      }

      outputFileU << "\n";

      //uSlice.clear();
    }      
    outputFileU.close();

    std::string filenameV = dir + "final/v.dat";
    std::ofstream outputFileV(filenameV);

    for(int i = 0; i < nx; ++i){
      int idy = (ny+1)/2;      
      if (exm == "tgv" || exm == "stgv"){
        idy = ny / 2 + 1;
      }
      //std::vector<double> vSlice(ops.No, 0.0);
      for (int alpha = 0; alpha < ops.No; ++alpha) {
        vSlice[alpha] = v[i][idy][alpha];
      }


      outputFileV.precision(20);  
      outputFileV << i * dx << "\t" << ops.mean(vSlice) * conversionVelocity << "\t" << ops.std(vSlice) * conversionVelocity;

      if (exm == "tgv" || exm == "stgv"){
        double x = i * dx;
        double y = idy * dy;
        double k2 = dx*dx + dy*dy;
        double damp = std::exp(-k2*physViscosity*iter);
        outputFileU << "\t" << u0 * std::sin(x) * std::cos(y) * damp * conversionVelocity;
      }

      outputFileV << "\n";
      // vSlice.clear();
    }      
    outputFileV.close();

  }
    
  void totalKineticEnergy(std::vector<double>&tke, double&tkeAna, int t)
    {
      std::vector<double> u2(ops.No, 0.0);
      std::vector<double> v2(ops.No, 0.0);
      std::vector<double> uSlice(ops.No, 0.0);
      std::vector<double> vSlice(ops.No, 0.0);
      std::vector<double> tkeChaos(ops.No, 0.0);
      
      // std::vector<double> u2Ran(ops.total_nq, 0.0);
      // std::vector<double> v2Ran(ops.total_nq, 0.0);
      // std::vector<double> tkeRan(ops.total_nq, 0.0);
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; ++j){
                for (int alpha = 0; alpha < ops.No; ++alpha) {
                    uSlice[alpha] = u[i][j][alpha];
                    vSlice[alpha] = v[i][j][alpha];
                }

                ops.chaos_product(uSlice, uSlice, u2);
                ops.chaos_product(vSlice, vSlice, v2);

                // ops.chaos2ran(u2, u2Ran);
                // ops.chaos2ran(v2, v2Ran);

                // for (int sample = 0; sample < ops.total_nq; ++sample) {
                //   tkeRan[sample] = ((u2Ran[sample] + v2Ran[sample]) *  2 / (nx*ny*u0*u0));
                // }

                // std::cout << "tkeRan: " << tkeRan[0] << std::endl;

                // ops.ran2chaos(tkeRan, tkeChaos);

                for (int alpha = 0; alpha < ops.No; ++alpha) {
                    tke[alpha] += ((u2[alpha] + v2[alpha]) *  2 / (nx*ny*u0*u0));
                    //tke[alpha] += tkeChaos[alpha];
                }
                
                double x = i * dx;
                double y = j * dy;
                double k2 = dx*dx + dy*dy;
                double damp = std::exp(-k2*physViscosity*t);
                double uAna = -u0 * std::cos(x) * std::sin(y) * damp;
                double vAna =  u0 * std::sin(x) * std::cos(y) * damp;
                tkeAna += ((uAna * uAna + vAna * vAna)*2/(nx*ny*u0*u0));

            }
        }
        // u2.clear();
        // v2.clear();
        // uSlice.clear();
        // vSlice.clear();
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
              //uSlice[alpha] = u[i][j][alpha];
              //vSlice[alpha] = v[i][j][alpha];
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
