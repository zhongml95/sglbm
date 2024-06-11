#include "tgv2d.h"

int main( int argc, char* argv[] )
{
  Parameters params;
  bool uq = false;

  // Call readParameters to populate the params instance
  readParameters("./parameters.dat", params);
  double physViscosity = params.physVelocity * params.L / params.Re;
  params.tau = 3 * physViscosity + 0.5;

  std::string dir    = "./data/nx" + std::to_string(params.resolution) + "/";
  std::string dirAna = "./data/nx" + std::to_string(params.resolution) + "/final/";

  std::string command;
  int a;
  command = "rm -rf " + dir;
  a = std::system(command.c_str());    
  command = "mkdir -p " + dir;
  a = std::system(command.c_str());
  command = "mkdir -p " + dirAna;
  a = std::system(command.c_str());
  
  std::cout << dir << std::endl;
  std::cout << "finish mkdir" << std::endl;
  double start = omp_get_wtime();

  simulateTGV2D(params, dir, 0, uq);

  double end = omp_get_wtime();
  std::cout << "total MCS time used: " <<  end - start << std::endl;

  return 0;
}