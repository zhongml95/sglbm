#include <cmath>
#include <random>
#include "../../src/lbm.h"

double calc_tke_error(lbm lbm, int count) {
  
  double tke = 0.0;
  double tkeAna = 0.0;
  lbm.totalKineticEnergy(tke, tkeAna, count);
  return std::abs((tke - tkeAna) / tkeAna);
}

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
  int resolution = 64;
  double L = 2.0 * M_PI;
  double lx = 2.0 * M_PI;
  double ly = 2.0 * M_PI;
  double dx = L / resolution;
  double dy = L / resolution;

  double Re = 15;
  double physVelocity = 0.01;
  double physViscosity = physVelocity*L/Re;
  double tau = 3 * physViscosity + 0.5;


  std::vector<std::vector<int>> material(resolution+1, std::vector<int>(resolution+1, 1));
    
  std::string dir    = "./data/nx" + std::to_string(resolution) + "/";
  std::string dirAna = "./data/nx" + std::to_string(resolution) + "/final/";

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

  lbm lbm(dir, "tgv");
  lbm.setGeometry(L,resolution,lx,ly,material);
  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;
  
  double start_mc = omp_get_wtime();
  
  lbm.setFluid(physVelocity, physViscosity, tau);
  lbm.initialize();
    
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

  std::cout << "total MCS time used: " <<  end - start_0 << std::endl;

  return 0;
}
