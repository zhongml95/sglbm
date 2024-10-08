#include <cmath>
#include "../../src/sglbm.h"

double calc_tke_error(sglbm sglbm, int count) {
  
  std::vector<double> tke(sglbm.ops.No,0.0);
  double tkeAna = 0.0;
  sglbm.totalKineticEnergy(tke, tkeAna, count);
  return std::abs((sglbm.ops.mean(tke)-tkeAna) / tkeAna);
}

int main( int argc, char* argv[] )
{
    Parameters params;
    
    // Call readParameters to populate the params instance
    readParameters("./parameters.dat", params);

    double dx = 2 * M_PI * params.L / params.resolution;
    double dy = 2 * M_PI * params.L / params.resolution;

    double dt = dx / (params.physVelocity / (params.Ma / std::sqrt(3)));

    double conversionViscosity = dx * dx  / dt;

    double physViscosity = params.physVelocity * params.L / params.Re;

    double tau = physViscosity / conversionViscosity * 3 + 0.5;
     

    std::vector<std::vector<int>> material(params.resolution+1, std::vector<int>(params.resolution+1, 1));
    
    //std::string dir = "./data/tgv/t5/ViscosityNr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(resolution) + "/";
    //std::string dirAna = "./data/tgv/t5/ViscosityNr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(resolution) + "/final/";
    
    std::string dir = "./data/tgv/t5/ReNr" + std::to_string(params.order) + "Nq" + std::to_string(params.nq) + "N" + std::to_string(params.resolution) + "/";
    std::string dirAna = "./data/tgv/t5/ReNr" + std::to_string(params.order) + "Nq" + std::to_string(params.nq) + "N" + std::to_string(params.resolution) + "/final/";

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


    sglbm sglbm(dir, "stgv", params);
    sglbm.setGeometry(params.L, params.resolution, params.lx, params.ly, material);
    sglbm.setFluid(params.physVelocity, physViscosity, tau);
    sglbm.initialize();
    //sglbm.iteration();


    std::cout << "start iteration" << std::endl;
    double maxPhysT = 30.0;
    int maxIter = maxPhysT / dt;
    std::cout << "maxIter: " << maxIter << std::endl;
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
    for (int i = 0; i < maxIter; ++i) {
      sglbm.collision();
      sglbm.streaming();      
      // #pragma omp single
      // {
      //     end = omp_get_wtime();
      sglbm.reconstruction();    
      
      //     std::cout << "reconstruction CPI time used: " << end - start << "s" << std::endl;   
      //     start = end;
      // }
#pragma omp single
      {
        if (i % 100 == 0) {
          //c_end = std::clock();
          end = omp_get_wtime();
          //double time_elapsed_s = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
          //std::cout << "iter: " << i << " " << "CPI time used: " << time_elapsed_s << "ms" << std::endl;
          
          err = calc_tke_error(sglbm, count);

          std::cout << "iter: " << i << " " << "CPI time used: " << end - start << "s" << "\t" << "TKE error " << err << std::endl;
          // sglbm.output(dir, i, end - start);
          //c_start = c_end;
          start = end;
          t = 0.0;
        }
      }
      count = i;
    }
    end = omp_get_wtime();
    sglbm.output(dir, maxIter, end - start_0);

    err = calc_tke_error(sglbm, count);

    std::cout << "total CPI time used: " << end - start_0 << "s" << "\t" << "TKE error " << err << std::endl;


    return 0;
}
