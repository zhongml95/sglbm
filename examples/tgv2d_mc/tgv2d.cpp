#include "../tgv2d_lbm/tgv2d.h"
#include <random>

int main( int argc, char* argv[] )
{
  Parameters params;
  bool uq = true;

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

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(params.parameter1[0], params.parameter2[0]);


  for ( int i = 0; i < params.nq; i++ ) {
    params.tau = 3 * (physViscosity * dis(gen)) + 0.5;
    
    simulateTGV2D(params, dir, i, uq);
  }

  return 0;
}