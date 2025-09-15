// postprocessing_sglbm.h
#ifndef POSTPROCESSING_SGLBM_H
#define POSTPROCESSING_SGLBM_H

#include "sglbm.h"

template<typename T>
void velocity_central(std::string dir, sglbm<T>& sglbm, UnitConverter<T>& uc, int idx, int idy) {

  // velocity in x direction at the center verticle line of the domain
  std::string filenameU = dir + "final/u.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < sglbm.ny; ++j){
      const size_t id = sglbm.idx(idx,j);
      outputFileU.precision(20);  
      outputFileU << j * sglbm.dx << "\t" << sglbm.ops->mean(sglbm.u[id]) * uc.getConversionVelocity() << "\t" << sglbm.ops->std(sglbm.u[id]) * uc.getConversionVelocity();

      outputFileU << "\n";
    }      
    outputFileU.close();


  // velocity in y direction at the center horizontal line of the domain
  std::string filenameV = dir + "final/v.dat";
  std::ofstream outputFileV(filenameV);
  for(int i = 0; i < sglbm.nx; ++i){
      const size_t id = sglbm.idx(i,idy);
      outputFileV.precision(20);  
      outputFileV << i * sglbm.dx << "\t" << sglbm.ops->mean(sglbm.v[id]) * uc.getConversionVelocity() << "\t" << sglbm.ops->std(sglbm.v[id]) * uc.getConversionVelocity();

      outputFileV << "\n";
    }      
    outputFileV.close();
}

template<typename T>
void velocity_all(std::string dir, sglbm<T>& sglbm, UnitConverter<T>& uc) {

  // mean velocity in x direction of the domain
  std::string filenameU = dir + "final/u_mean.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < sglbm.ny; ++j){
    for(int i = 0; i < sglbm.nx; ++i){
      const size_t id = sglbm.idx(i,j);
      outputFileU.precision(20);  
      outputFileU << sglbm.ops->mean(sglbm.u[id]) * uc.getConversionVelocity() << "\t" ;
    }      
    outputFileU << "\n";
  }
  outputFileU.close();

  // mean velocity in y direction of the domain
  std::string filenameV = dir + "final/v_mean.dat";
  std::ofstream outputFileV(filenameV);
  for(int j = 0; j < sglbm.ny; ++j) {
    for(int i = 0; i < sglbm.nx; ++i) {
      const size_t id = sglbm.idx(i,j);
      outputFileV.precision(20);  
      outputFileV << sglbm.ops->mean(sglbm.v[id]) * uc.getConversionVelocity() << "\t";
    }      
    outputFileV << "\n";
  }
  outputFileV.close();

  // std velocity in x direction of the domain
  std::string filename_u_std = dir + "final/u_std.dat";
  std::ofstream outputFileUStd(filename_u_std);
  for(int j = 0; j < sglbm.ny; ++j){
    for(int i = 0; i < sglbm.nx; ++i){
      const size_t id = sglbm.idx(i,j);
      outputFileUStd.precision(20);
      outputFileUStd << sglbm.ops->std(sglbm.u[id]) * uc.getConversionVelocity() << "\t" ;
    }      
    outputFileUStd << "\n";
  }
  outputFileUStd.close();

  // std velocity in y direction of the domain
  std::string filename_v_std = dir + "final/v_std.dat";
  std::ofstream outputFileVStd(filename_v_std);
  for(int j = 0; j < sglbm.ny; ++j) {
    for(int i = 0; i < sglbm.nx; ++i) {
      const size_t id = sglbm.idx(i,j);
      outputFileVStd.precision(20);  
      outputFileVStd << sglbm.ops->std(sglbm.v[id]) * uc.getConversionVelocity() << "\t";
    }      
    outputFileVStd << "\n";
  }
  outputFileVStd.close();
}

template<typename T>
void totalKineticEnergy(const sglbm<T>& sglbm, UnitConverter<T>& uc, std::vector<double>& tke, double& tkeAna, int t)
{
  double u0 = uc.getPhysVelocity();
  std::vector<double> u2Chaos(sglbm.No, 0.0);
  std::vector<double> v2Chaos(sglbm.No, 0.0);
  std::vector<double> tkeChaos(sglbm.No, 0.0);
  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      const size_t id = sglbm.idx(i,j);
      sglbm.ops->chaosProduct(sglbm.u[id], sglbm.u[id], u2Chaos);
      sglbm.ops->chaosProduct(sglbm.v[id], sglbm.v[id], v2Chaos);
      
      for (int alpha = 0; alpha < sglbm.No; ++alpha) {
        tke[alpha] += ((u2Chaos[alpha] + v2Chaos[alpha]) *  0.5 / (sglbm.nx*sglbm.ny*u0*u0));
      }
               
      double x = i * uc.getDx();
      double y = j * uc.getDx();
      double k2 = 2.0 * uc.getDx() * uc.getDx();
      double damp = std::exp(-k2 * uc.getPhysViscosity() * t);
      double uAna = -u0 * std::cos(x) * std::sin(y) * damp;
      double vAna =  u0 * std::sin(x) * std::cos(y) * damp;
      tkeAna += ((uAna * uAna + vAna * vAna) * 0.5 /(sglbm.nx * sglbm.ny * u0 * u0));

    }
  }
}

template<typename T>
void outputTKE(std::string dir, sglbm<T>& sglbm, UnitConverter<T>& uc, int t, double total_computational_time)
{
  std::string filenameTKE = dir + "final/tke.dat";
  std::ofstream outputFileTKE(filenameTKE);
  double tkeAna = 0.0;

  std::vector<double> tke(sglbm.No, 0.0);
  totalKineticEnergy(sglbm, uc, tke, tkeAna, t);

  outputFileTKE.precision(20);  
  std::cout << "tke: " << sglbm.ops->mean(tke) << " " << sglbm.ops->std(tke) << " " << tkeAna << " " << total_computational_time << std::endl;
  outputFileTKE << sglbm.ops->mean(tke) << "\t" << sglbm.ops->std(tke) << "\t" << tkeAna << "\t" << total_computational_time;
  outputFileTKE.close();
}


#endif // POSTPROCESSING_SGLBM_H