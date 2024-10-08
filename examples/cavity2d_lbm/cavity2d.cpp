#include "cavity2d.h"

int main( int argc, char* argv[] )
{
  Parameters params;
  bool uq = false;

  // Call readParameters to populate the params instance
  readParameters("./parameters.dat", params);

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

  simulateCavity2D(params, dir, 0, uq);

  return 0;
}