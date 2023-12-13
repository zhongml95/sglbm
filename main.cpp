#include <iostream>
#include <cmath>
//#include "src/polynomial.h"
#include "src/sglbm.h"

double calc_tke_error(sglbm sglbm, int count) {
  
  std::vector<double> tke(sglbm.op.order+1,0.0);
  double tkeAna = 0.0;
  sglbm.totalKineticEnergy(tke, tkeAna, count);
  return std::abs((sglbm.op.mean(tke)-tkeAna) / tkeAna);
}

int main( int argc, char* argv[] )
{
    unsigned int order = 4;
    unsigned int nq = 2*order+1;
    int polynomialType = 0;
    int resolution = 32;
    double parameter1 = 0.8;//0.8
    double parameter2 = 1.2;//1.2
    double L = 2.0 * M_PI;
    double lx = 2.0 * M_PI;
    double ly = 2.0 * M_PI;
    double dx = L / resolution;
    double dy = L / resolution;

    double Re = 15;
    double physVelocity = 0.01;
    double nu = physVelocity*L/Re;
    double tau = 3 * nu + 0.5;
    std::vector<std::vector<int>> material(resolution+1, std::vector<int>(resolution+1, 1));
    
    std::string dir = "./data/tgv/t5/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(resolution) + "/";
    std::string dirAna = "./data/tgv/t5/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(resolution) + "/final/";

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


    sglbm sglbm(dir, "tgv", nq, order, parameter1, parameter2, polynomialType);
    sglbm.setGeometry(L,resolution,lx,ly,material);
    sglbm.setFluid(physVelocity,nu,tau);
    sglbm.initialize();
    //sglbm.iteration();


    std::cout << "start iteration" << std::endl;
    double td = 1.0 / (sglbm.physViscosity * (dx * dx + dy * dy));
    std::cout << "td: " << td << std::endl;
    int count = 0;
    //std::clock_t c_start = std::clock();
    double start = omp_get_wtime();
    double start_0 = start;
    //std::clock_t c_end = std::clock();
    double end = omp_get_wtime();
    double err = 0.;

    double t = 0.0, t0, t1;

    sglbm.output(dir, 0, 0);


    size_t cores = omp_get_num_procs();
    omp_set_dynamic(0);
    omp_set_num_threads(cores);
    std::cout << "num Threads: " << cores << std::endl;
#pragma omp parallel 
    for (int i = 1; i < int(td * 0.5); ++i) {
    //for (int i = 1; i < 3; ++i) {
      sglbm.collision();
      //sglbm.boundary();
      sglbm.streaming();
      sglbm.reconstruction(); 

/*#pragma omp single
      {
        count = i;
        if (i % 1000 == 0) {
          //c_end = std::clock();
          end = omp_get_wtime();
          //double time_elapsed_s = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
          //std::cout << "iter: " << i << " " << "CPI time used: " << time_elapsed_s << "ms" << std::endl;
          
          err = calc_tke_error(sglbm, count);

          std::cout << "iter: " << i << " " << "CPI time used: " << end - start << "s" << "\t" << "TKE error " << err << std::endl;
          sglbm.output(dir, i, end - start);
          //c_start = c_end;
          start = end;
          t = 0.0;
        }
      }*/
    }
    count = int(td * 0.5) - 1;
    end = omp_get_wtime();
    sglbm.output(dir, count, end - start_0);

    err = calc_tke_error(sglbm, count);

    std::cout << "total CPI time used: " << end - start_0 << "s" << "\t" << "TKE error " << err << std::endl;


    return 0;
}
