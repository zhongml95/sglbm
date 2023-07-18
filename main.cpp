#include <iostream>
#include <cmath>
//#include "src/polynomial.h"
#include "src/sglbm.h"


int main( int argc, char* argv[] )
{
    unsigned int order = 15;
    unsigned int nq = 31;
    int resolution = 32;
    double parameter1 = 0.9;
    double parameter2 = 1.1;
    double L = 2.0 * M_PI;
    double lx = 2.0 * M_PI;
    double ly = 2.0 * M_PI;
    double dx = L / resolution;
    double dy = L / resolution;

    double tau = 0.5125663706143592;
    double Re = 15;
    double physVelocity = 0.01;
    double nu = physVelocity*2*M_PI/Re;
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


    sglbm sglbm(dir, "tgv", nq, order, parameter1,parameter2);
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
    //std::clock_t c_end = std::clock();
    double end = omp_get_wtime();

    double t = 0.0, t0, t1;

    sglbm.output(dir, 0);


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
#pragma omp single
      {
        count = i;
        if (i % 1000 == 0) {
          //c_end = std::clock();
          end = omp_get_wtime();
          //double time_elapsed_s = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
          //std::cout << "iter: " << i << " " << "CPI time used: " << time_elapsed_s << "ms" << std::endl;
          std::vector<double> tke(order+1,0.0);
          double tkeAna = 0.0;
          sglbm.totalKineticEnergy(tke, tkeAna, count);
          double err = std::abs((sglbm.mean(tke)-tkeAna) / tkeAna);
          std::cout << "iter: " << i << " " << "CPI time used: " << end - start << "s" << "\t" << "TKE error " << err << std::endl;
          //sglbm.output(dir, i);
          //c_start = c_end;
          start = end;
          t = 0.0;
        }
      }
    }
    sglbm.output(dir, count);


    return 0;
}