#ifndef LBM_H
#define LBM_H

#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <sys/stat.h>
#include <iomanip>
#include "util.h"

class lbm{
public:
  int N;
  int nx;
  int ny;
  double L;
  double dx;
  double lx;
  double ly;
  double physVelocity;
  double Re;
  double physViscosity;
  double physDensity;
  double tau;
  double dt;
  double conversionViscosity;
  double conversionVelocity;
  double conversionDensity;
  double conversionMass;
  double conversionForce;
  double u0;
  double omega;
  std::string dir;

  double cs2 = 1.0 / 3.0;
  std::vector<int> cx = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
  std::vector<int> cy = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
  std::vector<double> w = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
  std::vector<int> kinv = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

  std::vector<std::vector<int>> material;
  std::vector<std::vector<std::vector<double>>> bouzidiQ;
  std::vector<std::vector<double>> u;
  std::vector<std::vector<double>> v;
  std::vector<std::vector<double>> rho;
  
  std::vector<std::vector<std::vector<double>>> f;
  std::vector<std::vector<std::vector<double>>> F;
  std::vector<std::vector<std::vector<double>>> feq;

  lbm(std::string _dir)
  {
    dir = _dir;
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
    omega = 1.0 / (3 * physViscosity / conversionViscosity + 0.5);

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
        u[i][j] = 0.0;
        v[i][j] = 0.0;
        rho[i][j] = 0.0;
        
        f[i][j].resize(9);
        F[i][j].resize(9);
        feq[i][j].resize(9);
        for (int k = 0; k < 9; ++k) {
          f[i][j][k] = 0.0;
          F[i][j][k] = 0.0;
          feq[i][j][k] = 0.0;
        }
      }
    }

  }

  void initializeDistributionFunction() {

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          feq[i][j][k] = equilibrium(rho[i][j], u[i][j], v[i][j], cx[k], cy[k], w[k]);
          f[i][j][k] = feq[i][j][k];
          F[i][j][k] = feq[i][j][k];
        }
      }
    }

  }

  void collision()
  {
    //std::cout << "collision start" << std::endl;
#pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < 9; ++k) {
          F[i][j][k] = f[i][j][k] - (f[i][j][k] - feq[i][j][k]) * omega;
        }
      }
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

#pragma omp for collapse(2)
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        rho[i][j] = 0;
        double rhou = 0;
        double rhov = 0;

        for (int k = 0; k < 9; ++k) {
          rho[i][j] += f[i][j][k];
        }

        if (material[i][j] == 1) {
          //std::cout << "check" << std::endl;
          for (int k = 0; k < 9; ++k) {
            rhou += f[i][j][k] * cx[k];
            rhov += f[i][j][k] * cy[k];
          }

          u[i][j] = rhou / rho[i][j];
          v[i][j] = rhov / rho[i][j];

        }
        else if (material[i][j] == 2)
        {
          u[i][j] = 0.0;
          v[i][j] = 0.0;
        }
        else if (material[i][j] == 3){
          u[i][j] = u0;
          v[i][j] = 0.0;
        }
        for (int k = 0; k < 9; ++k) {
          feq[i][j][k] = equilibrium(rho[i][j], u[i][j], v[i][j], cx[k], cy[k], w[k]); 
        }
      } //
    }

  }

    void output(std::string dir, int iter, bool uq)
    {
      std::string filename;
      if (!uq) {
        filename = dir + "final.dat";
      }
      else {
        filename = dir + std::to_string(iter) + ".dat";
      }
      std::ofstream outputFile(filename);
      if (!outputFile) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
      }

      outputFile << "variables = \"X\", \"Y\", \"Rho\", \"Ux\", \"Vy\", \"GEOMETRY\"\n";
      outputFile << "ZONE I = " << nx << ", J = " << ny << ", F = POINT\n";

      for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
          double x = i * dx;
          double y = j * dx;
            
          outputFile << x << "\t" << y << "\t" << rho[i][j] << "\t" << u[i][j] * conversionVelocity << "\t" << v[i][j] * conversionVelocity << "\t" << material[i][j] <<  "\n";
        }
      }
      outputFile.close();

  }

  void boundary()
  {
      #pragma omp for collapse(2)
      for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
          if (material[i][j] == 2){
            // bounce back
            for (int k = 0; k < 9; ++k){
              int new_i = i+cx[k];
              int new_j = j+cy[k];
              if ((new_i!=-1) && (new_i!=nx) && (new_j!=0) && (new_j!=ny)){
                if (material[new_i][new_j] == 1){
                  F[i][j][k] = F[new_i][new_j][kinv[k]];
                }
              }
            }
          }
          else if (material[i][j] == 3){
            // moving wall
            double rhoWall = rho[i][j-1];
            for (int k = 0; k < 9; ++k) {
                F[i][j][k] = equilibrium(rhoWall, u0, 0.0, cx[k], cy[k], w[k]) + F[i][j-1][k] - feq[i][j-1][k];
            }
          }
        }
      }
    }

};



#endif // LBM_H
