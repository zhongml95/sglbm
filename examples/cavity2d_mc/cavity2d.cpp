#include "../cavity2d_lbm/cavity2d.h"
#include <random>

int main( int argc, char* argv[] )
{
  Parameters params;
  bool uq = true;

  // Call readParameters to populate the params instance
  readParameters("./parameters.dat", params);
  double physVelocity = params.physVelocity;

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
    params.physVelocity = physVelocity * dis(gen);
    
    simulateCavity2D(params, dir, i, uq);
  }

  return 0;
}