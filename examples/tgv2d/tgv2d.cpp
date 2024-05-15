#include <cmath>
#include "../../src/sglbm.h"
#include "../../src/postprocessing_sglbm.h"

// #define stochastic_Re
#define stochastic_viscosity

double calc_tke_error(sglbm sglbm, int count) {
  
  std::vector<double> tke(sglbm.ops.order+1,0.0);
  double tkeAna = 0.0;
  totalKineticEnergy(sglbm, tke, tkeAna, count);
  return std::abs((sglbm.ops.mean(tke)-tkeAna) / tkeAna);
}

void setGeometry(sglbm& sglbm, Parameters params) {
  std::cout << "start setting geometry" << std::endl;

  sglbm.nx = sglbm.N;
  sglbm.ny = sglbm.N;

  sglbm.material = std::vector<std::vector<int>>(sglbm.nx, std::vector<int>(sglbm.ny, 1));

  std::cout << "nx: " << sglbm.nx << ", ny: " << sglbm.ny << std::endl;

  std::cout << "finish setting geometry" << std::endl;

}

void initialize(sglbm& sglbm) {

  sglbm.prepareLattice();


  std::vector<double> omegaRan(sglbm.ops.total_nq, 0.0);


  #if defined(stochastic_Re)
    std::vector<double> ReChaos(sglbm.ops.No, 0.0);
    std::vector<double> ReRan(sglbm.ops.total_nq, 0.0);
    std::vector<double> chaos(2, 0.0);
    sglbm.ops.convert2affinePCE(sglbm.Re * sglbm.ops.parameter1[0], sglbm.Re * sglbm.ops.parameter2[0], sglbm.ops.polynomial_types[0], chaos);
    ReChaos[0] = chaos[0];
    ReChaos[1] = chaos[1];
    sglbm.ops.chaos2ran(ReChaos, ReRan);

    for (int sample = 0; sample < sglbm.ops.total_nq; ++sample) {
      omegaRan[sample] = 1.0 / ( 3.0 * (sglbm.physVelocity * sglbm.L / ReRan[sample]) / sglbm.conversionViscosity + 0.5 );
    }
  #elif defined(stochastic_viscosity)
    std::vector<double> physViscosityChaos(sglbm.ops.No, 0.0);
    std::vector<double> physViscosityRan(sglbm.ops.total_nq, 0.0);
    std::vector<double> chaos(2, 0.0);

    sglbm.ops.convert2affinePCE(sglbm.physViscosity * sglbm.ops.parameter1[0], sglbm.physViscosity * sglbm.ops.parameter2[0], sglbm.ops.polynomial_types[0], chaos);
    physViscosityChaos[0] = chaos[0];
    physViscosityChaos[1] = chaos[1];
    sglbm.ops.chaos2ran(physViscosityChaos, physViscosityRan);
          
    for (int sample = 0; sample < sglbm.ops.total_nq; ++sample) {
      omegaRan[sample] = 1.0 / ( 3.0 * physViscosityRan[sample] / sglbm.conversionViscosity + 0.5 );
    }
  #endif
  
  std::vector<double> omegaChaos(sglbm.ops.No, 0.0);
  sglbm.ops.ran2chaos(omegaRan, omegaChaos);
  sglbm.omegaChaos = omegaChaos;

  for (int i = 0; i < sglbm.ops.No; ++i) {
    std::cout << "omegaChaos[" << i << "]: " << sglbm.omegaChaos[i] << std::endl;
  }

  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      double x = i * sglbm.dx;
      double y = j * sglbm.dx;

      std::vector<double> uChaos(sglbm.ops.No, 0.0);
      std::vector<double> vChaos(sglbm.ops.No, 0.0);
      std::vector<double> rChaos(sglbm.ops.No, 0.0);

      rChaos[0] = 1.0 - 1.5 * sglbm.u0 * sglbm.u0 * std::cos(x + y) * std::cos(x - y);
      uChaos[0] = -sglbm.u0 * std::cos(x) * std::sin(y);
      vChaos[0] =  sglbm.u0 * std::sin(x) * std::cos(y);

      sglbm.rho[i][j] = rChaos;
      sglbm.u[i][j] = uChaos;
      sglbm.v[i][j] = vChaos;      
    }
  }

  sglbm.initializeDistributionFunction();

  std::cout << "finish initializing" << std::endl;
}


int main( int argc, char* argv[] )
{
  Parameters params;
    
  // Call readParameters to populate the params instance
  readParameters("./parameters.dat", params);
  double physViscosity = params.physVelocity * params.L / params.Re;
  params.tau = 3 * physViscosity + 0.5;

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


  sglbm sglbm(dir, params);
  sglbm.UnitConverterFromResolutionAndRelaxationTime(params);

  setGeometry(sglbm, params);

  initialize(sglbm);

  std::cout << "start iteration" << std::endl;
  double td = 1.0 / (sglbm.physViscosity * (sglbm.dx * sglbm.dx * 2.0));
  std::cout << "td: " << td << std::endl;
  int count = 0;
  //std::clock_t c_start = std::clock();
  double start = omp_get_wtime();
  double start_0 = start;
  //std::clock_t c_end = std::clock();
  double end = omp_get_wtime();
  double err = 0.;

  double t = 0.0, t0, t1;

  sglbm.output(dir, 0);


  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;
#pragma omp parallel 
  for (int i = 1; i < int(td * 0.5); ++i) {
    sglbm.collision();
    sglbm.streaming();
    sglbm.reconstruction(); 

#pragma omp single
    {
      if (i % 1000 == 0) {
        //c_end = std::clock();
        end = omp_get_wtime();
        //double time_elapsed_s = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
        //std::cout << "iter: " << i << " " << "CPI time used: " << time_elapsed_s << "ms" << std::endl;
          
        err = calc_tke_error(sglbm, count);

        std::cout << "iter: " << i << " " << "CPI time used: " << end - start << "s" << "\t" << "TKE error " << err << std::endl;
        sglbm.output(dir, i);

        //c_start = c_end;
        start = end;
        t = 0.0;
      }
    }
    count = i;
  }
  count = int(td * 0.5) - 1;
  end = omp_get_wtime();
  sglbm.output(dir, count);

  err = calc_tke_error(sglbm, count);

  std::cout << "total CPI time used: " << end - start_0 << "s" << "\t" << "TKE error " << err << std::endl;

  velocity_central(dir, sglbm, (sglbm.nx/2)+1, (sglbm.ny/2)+1);
  velocity_all(dir, sglbm);
  outputTKE(dir, sglbm, count, end - start_0);



  return 0;
}
