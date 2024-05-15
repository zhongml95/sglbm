#include "sglbm.h"

void velocity_central(std::string dir, sglbm& sglbm, int idx, int idy) {

  // velocity in x direction at the center verticle line of the domain
  std::string filenameU = dir + "final/u.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < sglbm.ny; ++j){

      outputFileU.precision(20);  
      outputFileU << j * sglbm.dx << "\t" << sglbm.ops.mean(sglbm.u[idx][j]) * sglbm.conversionVelocity << "\t" << sglbm.ops.std(sglbm.u[idx][j]) * sglbm.conversionVelocity;

      outputFileU << "\n";
    }      
    outputFileU.close();


  // velocity in y direction at the center horizontal line of the domain
  std::string filenameV = dir + "final/v.dat";
  std::ofstream outputFileV(filenameV);
  for(int i = 0; i < sglbm.nx; ++i){

      outputFileV.precision(20);  
      outputFileV << i * sglbm.dx << "\t" << sglbm.ops.mean(sglbm.v[i][idy]) * sglbm.conversionVelocity << "\t" << sglbm.ops.std(sglbm.v[i][idy]) * sglbm.conversionVelocity;

      outputFileV << "\n";
    }      
    outputFileV.close();
}

void velocity_all(std::string dir, sglbm& sglbm) {

  // mean velocity in x direction of the domain
  std::string filenameU = dir + "final/u_mean.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < sglbm.ny; ++j){
    for(int i = 0; i < sglbm.nx; ++i){
      outputFileU.precision(20);  
      outputFileU << sglbm.ops.mean(sglbm.u[i][j]) * sglbm.conversionVelocity << "\t" ;
    }      
    outputFileU << "\n";
  }
  outputFileU.close();

  // mean velocity in y direction of the domain
  std::string filenameV = dir + "final/v_mean.dat";
  std::ofstream outputFileV(filenameV);
  for(int j = 0; j < sglbm.ny; ++j) {
    for(int i = 0; i < sglbm.nx; ++i) {
      outputFileV.precision(20);  
      outputFileV << sglbm.ops.mean(sglbm.v[i][j]) * sglbm.conversionVelocity << "\t";
    }      
    outputFileV << "\n";
  }
  outputFileV.close();

  // std velocity in x direction of the domain
  std::string filenameU = dir + "final/u_std.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < sglbm.ny; ++j){
    for(int i = 0; i < sglbm.nx; ++i){
      outputFileU.precision(20);  
      outputFileU << sglbm.ops.std(sglbm.u[i][j]) * sglbm.conversionVelocity << "\t" ;
    }      
    outputFileU << "\n";
  }
  outputFileU.close();

  // std velocity in y direction of the domain
  std::string filenameV = dir + "final/v_std.dat";
  std::ofstream outputFileV(filenameV);
  for(int j = 0; j < sglbm.ny; ++j) {
    for(int i = 0; i < sglbm.nx; ++i) {
      outputFileV.precision(20);  
      outputFileV << sglbm.ops.std(sglbm.v[i][j]) * sglbm.conversionVelocity << "\t";
    }      
    outputFileV << "\n";
  }
  outputFileV.close();
}


void totalKineticEnergy(sglbm& sglbm, std::vector<double>&tke, double&tkeAna, int t)
{
  std::vector<double> u2Chaos(sglbm.ops.No, 0.0);
  std::vector<double> v2Chaos(sglbm.ops.No, 0.0);
  std::vector<double> tkeChaos(sglbm.ops.No, 0.0);
      
  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      sglbm.ops.chaos_product(sglbm.u[i][j], sglbm.u[i][j], u2Chaos);
      sglbm.ops.chaos_product(sglbm.v[i][j], sglbm.v[i][j], v2Chaos);
      
      for (int alpha = 0; alpha < sglbm.ops.No; ++alpha) {
        tke[alpha] += ((u2Chaos[alpha] + v2Chaos[alpha]) *  0.5 / (sglbm.nx*sglbm.ny*sglbm.u0*sglbm.u0));
      }
               
      double x = i * sglbm.dx;
      double y = j * sglbm.dx;
      double k2 = 2.0 * sglbm.dx * sglbm.dx;
      double damp = std::exp(-k2 * sglbm.physViscosity * t);
      double uAna = -sglbm.u0 * std::cos(x) * std::sin(y) * damp;
      double vAna =  sglbm.u0 * std::sin(x) * std::cos(y) * damp;
      tkeAna += ((uAna * uAna + vAna * vAna) * 0.5 /(sglbm.nx * sglbm.ny * sglbm.u0 * sglbm.u0));

    }
  }
}


void outputTKE(std::string dir, sglbm& sglbm, int t, double total_computational_time)
{
  std::string filenameTKE = dir + "final/tke.dat";
  std::ofstream outputFileTKE(filenameTKE);
  double tkeAna = 0.0;

  std::vector<double> tke(sglbm.ops.No, 0.0);
  totalKineticEnergy(sglbm, tke, tkeAna, t);

  outputFileTKE.precision(20);  
  std::cout << "tke: " << sglbm.ops.mean(tke) << " " << sglbm.ops.std(tke) << " " << tkeAna << " " << total_computational_time << std::endl;
  outputFileTKE << sglbm.ops.mean(tke) << "\t" << sglbm.ops.std(tke) << "\t" << tkeAna << "\t" << total_computational_time;
  outputFileTKE.close();
}
