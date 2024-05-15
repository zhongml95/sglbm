#include "../../src/sglbm.h"
#include "../../src/postprocessing_sglbm.h"

double calcError(sglbm sglbm, double&uNorm0) {
  double error = 0.0;
  double uNorm1 = 0.0;
  for (int i = 0; i < sglbm.nx; ++i){
    for (int j = 0; j < sglbm.ny; j++){
      std::vector<double> uSlice(sglbm.ops.order+1, 0.0);
      std::vector<double> vSlice(sglbm.ops.order+1, 0.0);
      for (int alpha = 0; alpha < sglbm.ops.order+1; ++alpha){
        uSlice[alpha] = sglbm.u[i][j][alpha];
        vSlice[alpha] = sglbm.v[i][j][alpha];
      }
      double uMean = 0.0;
      double vMean = 0.0;
      uMean = sglbm.ops.mean(uSlice);
      vMean = sglbm.ops.mean(vSlice);
      uNorm1 += (uMean * uMean + vMean * vMean);
    }
  }

  uNorm1 = std::sqrt(uNorm1);
  error = std::fabs(uNorm1 - uNorm0) / uNorm0;
  uNorm0 = uNorm1;
  return error;
}

void setGeometry(sglbm& sglbm, Parameters params) {
  std::cout << "start setting geometry" << std::endl;

  sglbm.nx = sglbm.N + 1;
  sglbm.ny = sglbm.N + 1;

  sglbm.material = std::vector<std::vector<int>>(sglbm.nx, std::vector<int>(sglbm.ny, 1));

  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      if (j == sglbm.ny-1) {
        sglbm.material[i][j] = 3;
      }
      
      if ((i == 0) || (i == sglbm.nx-1) || (j == 0)) {
        sglbm.material[i][j] = 2;
      }
    }
  }

  std::cout << "finish setting geometry" << std::endl;
    
}

void initialize(sglbm& sglbm) {

  sglbm.prepareLattice();

  std::vector<double> omegaChaos(sglbm.ops.No+1, 0.0);
  omegaChaos[0] = sglbm.omega0;
  sglbm.omegaChaos = omegaChaos;

  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      std::vector<double> uChaos(sglbm.ops.No+1, 0.0);
      std::vector<double> vChaos(sglbm.ops.No+1, 0.0);
      std::vector<double> rChaos(sglbm.ops.No+1, 0.0);

      rChaos[0] = 1.0;
      sglbm.rho[i][j] = rChaos;
      sglbm.u[i][j] = uChaos;
      sglbm.v[i][j] = vChaos;
    }
  }
  
  sglbm.initializeDistributionFunction();

  std::cout << "finish initializing" << std::endl;
}

void setBoundaryValue(sglbm& sglbm) {

  std::vector<double> chaos(2,0.0);
  sglbm.ops.convert2affinePCE(sglbm.ops.parameter1[0]*sglbm.u0, sglbm.ops.parameter2[0]*sglbm.u0, sglbm.ops.polynomial_types[0],chaos);

  for (int i = 0; i < sglbm.nx; ++i) {
    for (int j = 0; j < sglbm.ny; ++j) {
      if (sglbm.material[i][j] == 3) {
        sglbm.u[i][j][0] = chaos[0];
        sglbm.v[i][j][0] = chaos[1];
      }
    }
  }
}

int main( int argc, char* argv[] )
{
    Parameters params;
    
    // Call readParameters to populate the params instance
    readParameters("./parameters.dat", params);

    std::string dir = "./data/cavity2d/Nr" + std::to_string(params.order) + "Nq" + std::to_string(params.nq) + "N" + std::to_string(params.resolution) + "/";
    std::string dirAna = "./data/cavity2d/Nr" + std::to_string(params.order) + "Nq" + std::to_string(params.nq) + "N" + std::to_string(params.resolution) + "/final/";

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

    setBoundaryValue(sglbm);
    //sglbm.iteration();

    std::cout << "start iteration" << std::endl;
    int count = 1;
    //std::clock_t c_start = std::clock();
    double start = omp_get_wtime();
    //std::clock_t c_end = std::clock();
    double end = omp_get_wtime();

    double t = 0.0, t0, t1;

    sglbm.output(dir, 0);


    size_t cores = omp_get_num_procs();
    omp_set_dynamic(0);
    omp_set_num_threads(8);

  double error = 1.0;
  double uNorm = 1.0;


#pragma omp parallel 
      while(error > 0.0001) {
      sglbm.collision();
      sglbm.boundary();
      sglbm.streaming();
      sglbm.reconstruction(); // parallel for
#pragma omp single
      {
        count++;
        if (count % 1000 == 0){
          error = calcError(sglbm, uNorm);
          end = omp_get_wtime();
          std::cout << "iter: " << count << " " << "CPI time used: " << end - start << "s" << "\t" << "err: " << error << std::endl;
          start = end;
          t = 0.0;
          if (count % 10000 == 0) {
            sglbm.output(dir, count);
          }
        }
      }
    }
    sglbm.output(dir, count);
    velocity_central(dir, sglbm, (sglbm.nx+1)/2, (sglbm.ny+1)/2);
    velocity_all(dir, sglbm);
    return 0;
}

