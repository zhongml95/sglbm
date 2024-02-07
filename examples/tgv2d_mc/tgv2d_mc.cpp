#include <cmath>
//#include "src/polynomial.h"
#include "../../src/lbm.h"
#include "../../src/util.h"

double calc_tke_error(lbm lbm, int count) {
  
  double tke = 0.0;
  double tkeAna = 0.0;
  lbm.totalKineticEnergy(tke, tkeAna, count);
  return std::abs((tke - tkeAna) / tkeAna);
}

int main( int argc, char* argv[] )
{
    Parameters params;
    
    // Call readParameters to populate the params instance
    readParameters("./parameters.dat", params);

    double dx = params.L / params.resolution;
    double dy = params.L / params.resolution;

    double physViscosity = params.physVelocity * params.L / params.Re;
    double tau = 3 * physViscosity + 0.5;
    std::vector<std::vector<int>> material(params.resolution+1, std::vector<int>(params.resolution+1, 1));
    

  double parameter1 = 0.8 * physViscosity;//0.8
  double parameter2 = 1.2 * physViscosity;//1.2


  std::string dir    = "./data/tgv/t5/MC" + std::to_string(params.resolution) + "/";
  std::string dirAna = "./data/tgv/t5/MC" + std::to_string(params.resolution) + "/final/";

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

  int total_nq = 200;
  lbm lbm(dir, "tgv");
  lbm.setGeometry(params.L, params.resolution,params.lx,params.ly,material);
  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;
  
  // Seed the random number generator with the current time
  std::srand(static_cast<unsigned>(std::time(nullptr)));

  double start_mc = omp_get_wtime();
  for (int n = 0; n < total_nq; ++n) {  
    double random_physViscosity = parameter1 + static_cast<double>(std::rand()) / RAND_MAX * (parameter2 - parameter1);
    lbm.setFluid(params.physVelocity, random_physViscosity, tau);
    lbm.initialize();
    //sglbm.iteration();
    std::cout << "start iteration" << std::endl;
    double td = 1.0 / (lbm.physViscosity * (dx * dx + dy * dy));
    std::cout << "td: " << td << std::endl;
    int count = 0;
    //std::clock_t c_start = std::clock();
    double start = omp_get_wtime();
    double start_0 = start;
    //std::clock_t c_end = std::clock();
    double end = omp_get_wtime();
    double err = 0.;

    double t = 0.0, t0, t1;
    
  #pragma omp parallel 
    for (int i = 1; i < int(td * 0.5); ++i) {
      //for (int i = 1; i < 3; ++i) {
      lbm.collision();
      //sglbm.boundary();
      lbm.streaming();
      lbm.reconstruction(); 
  
    }

  end = omp_get_wtime();
  lbm.output(dir, int(td * 0.5) - 1, end - start_0);
      //save_velocity_field(dir, resolution, resolution, lbm.u, lbm.v, n);

      //end = omp_get_wtime();
      //err = calc_tke_error(lbm, int(td * 0.5) - 1);

      //std::cout << "total CPI time used: " << end - start_0 << "s" << "\t" << "TKE error " << err << std::endl;
  }

  double end_mc = omp_get_wtime();

  std::cout << "total MCS time used: " <<  end_mc - start_mc << std::endl;


    return 0;
}
