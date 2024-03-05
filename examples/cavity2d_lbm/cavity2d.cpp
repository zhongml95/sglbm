#include <cmath>
#include <random>
#include "../../src/lbm.h"
#include "../../src/util.h"


void save_velocity_field(std::string dir, int nx, int ny, std::vector<std::vector<double>> u, std::vector<std::vector<double>> v, int iter) {
    std::string filename_u = dir + "/u_" + std::to_string(iter) + ".dat";
    std::ofstream outputFile_u(filename_u);
    if (!outputFile_u) {
      std::cerr << "Error opening the file: " << filename_u << std::endl;
      return;
    }
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x = i ;
        double y = j ;
        
        outputFile_u << u[i][j] << "\t";
      }
      outputFile_u << "\n";
    }

    std::string filename_v = dir + "/v_" + std::to_string(iter) + ".dat";
    std::ofstream outputFile_v(filename_v);
    if (!outputFile_v) {
      std::cerr << "Error opening the file: " << filename_v << std::endl;
      return;
    }
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        double x = i ;
        double y = j ;
        
        outputFile_v << v[i][j] << "\t";
      }
      outputFile_v << "\n";
    }

  }

int main( int argc, char* argv[] )
{
    Parameters params;
    
    // Call readParameters to populate the params instance
    readParameters("./parameters.dat", params);

    double dx = params.L / params.resolution;
    double dy = params.L / params.resolution;
    int nx = int(params.lx / dx) + 1;
    int ny = int(params.ly / dy) + 1;

    double physViscosity = params.physVelocity * params.L / params.Re;
    double tau = 0.501536;

    std::vector<std::vector<int>> material(params.resolution+1, std::vector<int>(params.resolution+1, 1));

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            if (j == ny-1)
            {
                material[i][j] = 3;
            }
            if ((i == 0) || (i == nx-1) || (j == 0)){
                material[i][j] = 2;
            }
        }
    }
    
    
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

  lbm lbm(dir, "cavity2d");
  lbm.setGeometry(params.L,params.resolution,params.lx,params.ly,material);
  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;
  
  double start_mc = omp_get_wtime();
  
  lbm.setFluid(params.physVelocity, physViscosity, tau);
  lbm.initialize();
    
  std::cout << "start iteration" << std::endl;
  
  int count = 0;
  //std::clock_t c_start = std::clock();
  double start = omp_get_wtime();
  double start_0 = start;
  //std::clock_t c_end = std::clock();
  double end = omp_get_wtime();
  double error = 1.0;

    
#pragma omp parallel 
    //for (int i = 1; i < 100000; ++i) {
    while(count * lbm.dt < 1.0) {
      //for (int i = 1; i < 3; ++i) {
      lbm.collision();
      lbm.boundary();
      lbm.streaming();
      lbm.reconstruction(); 
      #pragma omp single
      {
        count++;
        if (count % 100 == 0){
          end = omp_get_wtime();
          std::cout << "iter: " << count << " " << "CPI time used: " << end - start << "s" << std::endl;
          start = end;
          if (count % 1000 == 0) {
            lbm.output(dir, count, 0);
          }
        }
      }
    }

  end = omp_get_wtime();
  lbm.output(dir, count, 0);

  std::cout << "total MCS time used: " <<  end - start_0 << std::endl;

  return 0;
}
