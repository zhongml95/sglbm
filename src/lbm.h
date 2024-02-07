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

class lbm{
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

  lbm(std::string _dir, std::string _exm)
  {
    dir = _dir;
    exm = _exm;
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

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        //TGV            
          //printf("i: %d, j: %d, alpha: %d, omp_get_thread_num: %d\n", i,j,alpha,omp_get_thread_num()); 
        if(exm == "tgv") {
            double x = i * dx;
            double y = j * dy;
              
            rho[i][j] = 1.0 - 1.5 * u0 * u0 * std::cos(x + y) * std::cos(x - y);
            u[i][j] = -u0 * std::cos(x) * std::sin(y);
            v[i][j] = u0 * std::sin(x) * std::cos(y);
          }
          else if (exm == "cavity2d"){
            rho[i][j] = 1.0;
            //omegaChaos[0] = 1.0 / (3 * physViscosity / conversionViscosity + 0.5);
          }
      }
    }
    
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

          F[i][j][k] = f[i][j][k] - (f[i][j][k] - feq[i][j][k]) / tau; 

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
        /*else if (material[i][j] == 2)
        {
          for (int alpha = 0; alpha < op.order + 1; ++alpha) {
            u[i][j][alpha] = 0.0;
            v[i][j][alpha] = 0.0;
          }
          for (int sample = 0; sample < op.nq; sample++) {
            uRan[sample] = 0.0;
            vRan[sample] = 0.0;
          }
        }
        else if (material[i][j] == 3){
          std::vector<double> uSlice(op.order + 1, 0.0);
          std::vector<double> vSlice(op.order + 1, 0.0);
          std::vector<double> chaos(2,0.0);
          op.convert2affinePCE(op.parameter1*u0, op.parameter2*u0,chaos);
          uSlice[0] = chaos[0];
          uSlice[1] = chaos[1];
          op.evaluatePCE(uSlice, uRan);
          std::vector<double> empty(op.nq, 0.0);
          vRan.swap(empty);

          for (int alpha = 0; alpha < op.order + 1; ++alpha) {
            if (alpha == 0 || alpha == 1)
            u[i][j][alpha] = 0.0;
            v[i][j][alpha] = 0.0;
          }
        }*/

        for (int k = 0; k < 9; ++k) {
          feq[i][j][k] = equilibrium(rho[i][j], u[i][j], v[i][j], cx[k], cy[k], w[k]); 
        }
      } //
    }

  }

  void output(std::string dir, int iter, double time_cost)
  {
    /*std::string filename = dir + std::to_string(iter) + ".dat";
    std::ofstream outputFile(filename);
    if (!outputFile) {
      std::cerr << "Error opening the file: " << filename << std::endl;
      return;
    }

    outputFile << "variables = \"X\", \"Y\", \"Rho\", \"Ux\", \"Vy\", \"GEOMETRY\"\n";
    //outputFile << "variables = \"X\", \"Y\", \"RhoMean\", \"uMean\", \"vMean\", \"RhoStd\", \"uStd\", \"vStd\", \"geometry\"\n";
    outputFile << "ZONE I = " << nx << ", J = " << ny << ", F = POINT\n";

    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x = i * dx;
        double y = j * dy;
        
        outputFile << x << "\t" << y << "\t" << rho[i][j] << "\t" << u[i][j] * conversionVelocity << "\t" << v[i][j] * conversionVelocity << "\t" << material[i][j] <<  "\n";
        //outputFile << x << "\t" << y << "\t" << rho[i][j] << "\t" << u[i][j] * conversionVelocity << "\t" << v[i][j] * conversionVelocity << "\t" << u[i][j] * conversionVelocity << "\t" << v[i][j] * conversionVelocity << "\t" << u[i][j] * conversionVelocity  << "\t" << material[i][j] <<  "\n";

      }
    }
    outputFile.close();*/

    if (exm == "tgv"){
      std::string filenameTKE = dir + "final/tke.dat";
      std::ofstream outputFileTKE(filenameTKE, std::ios::app);
      double tke = 0.0;
      double tkeAna = 0.0;
      totalKineticEnergy(tke, tkeAna, iter+1);

      outputFileTKE.precision(20);  
      outputFileTKE << tke << "\t" << tkeAna << "\t" << time_cost << std::endl;
      outputFileTKE.close();
    }

    std::string filenameU = dir + "final/u.dat";
    std::ofstream outputFileU(filenameU);

    for(int j = 0; j < ny; ++j){
      int idx = (nx+1)/2;
      if (exm == "tgv"){
        idx = nx / 2 + 1;
      }

      outputFileU.precision(20);  
      outputFileU << j * dy << "\t" << u[idx][j] * conversionVelocity;

      if (exm == "tgv"){
        double x = idx * dx;
        double y = j * dy;
        double k2 = dx*dx + dy*dy;
        double damp = std::exp(-k2*physViscosity*iter);
        outputFileU << "\t" << -u0 * std::cos(x) * std::sin(y) * damp * conversionVelocity;
      }

      outputFileU << "\n";

    }      
    outputFileU.close();

    std::string filenameV = dir + "final/v.dat";
    std::ofstream outputFileV(filenameV);

    for(int i = 0; i < nx; ++i){
      int idy = (ny+1)/2;      
      if (exm == "tgv"){
        idy = ny / 2 + 1;
      }

      outputFileV.precision(20);  
      outputFileV << i * dx << "\t" << v[i][idy] * conversionVelocity;

      if (exm == "tgv"){
        double x = i * dx;
        double y = idy * dy;
        double k2 = dx*dx + dy*dy;
        double damp = std::exp(-k2*physViscosity*iter);
        outputFileU << "\t" << u0 * std::sin(x) * std::cos(y) * damp * conversionVelocity;
      }

      outputFileV << "\n";
    }      
    outputFileV.close();

  }

    void totalKineticEnergy(double&tke, double&tkeAna, int t)
    {
      for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
          tke += ((u[i][j] * u[i][j] + v[i][j] * v[i][j]) *  2 / (nx*ny*u0*u0));
                
          double x = i * dx;
          double y = j * dy;
          double k2 = dx*dx + dy*dy;
          double damp = std::exp(-k2*physViscosity*t);
          double uAna = -u0 * std::cos(x) * std::sin(y) * damp;
          double vAna =  u0 * std::sin(x) * std::cos(y) * damp;
          tkeAna += ((uAna * uAna + vAna * vAna)*2/(nx*ny*u0*u0));

        }
      }
    }

    /*void boundary()
    {
      std::vector<double> rRan(op.nq, 0.0);
      std::vector<double> uRan(op.nq, 0.0);
      std::vector<double> vRan(op.nq, 0.0);
      std::vector<double> rSlice(op.order+1, 0.0);
      std::vector<double> uSlice(op.order+1, 0.0);
      std::vector<double> feqSlice(op.order+1, 0.0);
      std::vector<double> chaos(2,0.0);
      op.convert2affinePCE(op.parameter1*u0, op.parameter2*u0,chaos);
      uSlice[0] = chaos[0];
      uSlice[1] = chaos[1];
      op.evaluatePCE(uSlice, uRan);

      #pragma omp for collapse(2)
      for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
          if (material[i][j] == 2){
            for (int k = 0; k < 9; ++k){
              int new_i = i+cx[k];
              int new_j = j+cy[k];
              if ((new_i!=-1) && (new_i!=nx) && (new_j!=0) && (new_j!=ny)){
                if (material[new_i][new_j] == 1){
                  for (int alpha = 0; alpha < op.order+1; ++alpha){
                    F[i][j][k][alpha] = F[new_i][new_j][kinv[k]][alpha];
                  }
                }
              }
            }
          }
          else if (material[i][j] == 3){
            for (int alpha = 0; alpha < op.order + 1; ++alpha) {
              rSlice[alpha] = rho[i][j-1][alpha];
              //uSlice[alpha] = u[i][j][alpha];
              //vSlice[alpha] = v[i][j][alpha];
            }
            op.evaluatePCE(rSlice, rRan);

            for (int k = 0; k < 9; ++k) {
              std::vector<double> feqRan(op.nq, 0.0);
              for (int sample = 0; sample < op.nq; sample++) {
                feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
              }

              op.ran2chaos(feqRan, feqSlice);
              for (int alpha = 0; alpha < op.order + 1; ++alpha) {
                F[i][j][k][alpha] = feqSlice[alpha] + F[i][j-1][k][alpha] - feq[i][j-1][k][alpha];
              }
            }
          }
        }
      }
    }*/

};



#endif // LBM_H
