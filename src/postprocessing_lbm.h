#include "lbm.h"

void velocity_central(std::string dir, lbm& lbm, int idx, int idy) {

  // velocity in x direction at the center verticle line of the domain
  std::string filenameU = dir + "final/u.dat";
  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < lbm.ny; ++j){

      outputFileU.precision(20);  
      outputFileU << j * lbm.dx << "\t" << lbm.u[idx][j] * lbm.conversionVelocity << "\t" << lbm.u[idx][j] * lbm.conversionVelocity;

      outputFileU << "\n";
    }      
    outputFileU.close();


  // velocity in y direction at the center horizontal line of the domain
  std::string filenameV = dir + "final/v.dat";
  std::ofstream outputFileV(filenameV);
  for(int i = 0; i < lbm.nx; ++i){

      outputFileV.precision(20);  
      outputFileV << i * lbm.dx << "\t" << lbm.v[i][idy] * lbm.conversionVelocity << "\t" << lbm.v[i][idy] * lbm.conversionVelocity;

      outputFileV << "\n";
    }      
    outputFileV.close();
}


void velocity_All(std::string dir, lbm& lbm, int index, bool uq) {

  // velocity in x direction of the domain
  std::string filenameU, filenameV;
  if (uq) {
    filenameU = dir + "final/uAll_" + std::to_string(index) + ".dat";
    filenameV = dir + "final/vAll_" + std::to_string(index) + ".dat";
  }
  else {
    filenameU = dir + "final/uAll.dat";
    filenameV = dir + "final/vAll.dat";
  }

  std::ofstream outputFileU(filenameU);
  for(int j = 0; j < lbm.ny; ++j){
    for(int i = 0; i < lbm.nx; ++i){
      outputFileU.precision(20);  
      outputFileU << lbm.u[i][j] * lbm.conversionVelocity;

      outputFileU << "\t";
    }      
    outputFileU << "\n";
  }
  outputFileU.close();

  // velocity in y direction of the domain
  std::ofstream outputFileV(filenameV);
  for(int j = 0; j < lbm.ny; ++j) {
    for(int i = 0; i < lbm.nx; ++i) {
      outputFileV.precision(20);  
      outputFileV << lbm.v[i][j] * lbm.conversionVelocity;

      outputFileV << "\t";
    }      
    outputFileV << "\n";
  }
  outputFileV.close();
}

void totalKineticEnergy(lbm& lbm, double& tke, double& tkeAna, int t)
{
  for (int i = 0; i < lbm.nx; ++i){
    for (int j = 0; j < lbm.ny; ++j){
      tke += ((lbm.u[i][j] * lbm.u[i][j] + lbm.v[i][j] * lbm.v[i][j]) *  0.5 / (lbm.nx*lbm.ny*lbm.u0*lbm.u0));
                
      double x = i * lbm.dx;
      double y = j * lbm.dx;
      double k2 = lbm.dx*lbm.dx*2;
      double damp = std::exp(-k2*lbm.physViscosity*t);
      double uAna = -lbm.u0 * std::cos(x) * std::sin(y) * damp;
      double vAna =  lbm.u0 * std::sin(x) * std::cos(y) * damp;
      tkeAna += ((uAna * uAna + vAna * vAna)*0.5/(lbm.nx*lbm.ny*lbm.u0*lbm.u0));
    }
  }
}


void outputTKE(const std::string& dir, lbm& lbm, int t, double total_computational_time, bool uq) {
  
  std::string filenameTKE = dir + "final/tke.dat";
  std::ios_base::openmode mode = uq ? std::ios::app : std::ios::out;
    
  std::ofstream outputFileTKE(filenameTKE, mode);
  if (!outputFileTKE.is_open()) {
    std::cerr << "Failed to open file: " << filenameTKE << std::endl;
    return;
  }

  double tkeAna = 0.0;
  double tke = 0.0;
  totalKineticEnergy(lbm, tke, tkeAna, t);

  outputFileTKE.precision(20);  
  std::cout << "tke: " << tke << " " << tkeAna << " " << total_computational_time << std::endl;
  outputFileTKE << tke << "\t" << tkeAna << "\t" << total_computational_time << std::endl;

  outputFileTKE.close();

}
