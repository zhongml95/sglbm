#include <cmath>
#include <random>
#include "../../src/lbm.h"
#include "../../src/postprocessing_lbm.h"

double calcError(lbm lbm, double&uNorm0) {
  double error = 0.0;
  double uNorm1 = 0.0;
  for (int i = 0; i < lbm.nx; ++i){
    for (int j = 0; j < lbm.ny; j++){
      uNorm1 += (lbm.u[i][j] * lbm.u[i][j] + lbm.v[i][j] * lbm.v[i][j]);
    }
  }

  uNorm1 = std::sqrt(uNorm1);
  error = std::fabs(uNorm1 - uNorm0) / uNorm0;
  uNorm0 = uNorm1;
  return error;
}

void setGeometry(lbm& lbm, Parameters params) {
  std::cout << "start setting geometry" << std::endl;

  lbm.nx = lbm.N + 1;
  lbm.ny = lbm.N + 1;

  lbm.material = std::vector<std::vector<int>>(lbm.nx, std::vector<int>(lbm.ny, 1));

  for (int i = 0; i < lbm.nx; ++i) {
    for (int j = 0; j < lbm.ny; ++j) {
      if (j == lbm.ny-1) {
        lbm.material[i][j] = 3;
      }
      
      if ((i == 0) || (i == lbm.nx-1) || (j == 0)) {
        lbm.material[i][j] = 2;
      }
    }
  }


  std::cout << "nx: " << lbm.nx << ", ny: " << lbm.ny << std::endl;

  std::cout << "finish setting geometry" << std::endl;

}

void initialize(lbm& lbm) {

  lbm.prepareLattice();

  for (int i = 0; i < lbm.nx; ++i) {
    for (int j = 0; j < lbm.ny; ++j) {
      double x = i * lbm.dx;
      double y = j * lbm.dx;
              
      lbm.rho[i][j] = 1.0 - 1.5 * lbm.u0 * lbm.u0 * std::cos(x + y) * std::cos(x - y);
      lbm.u[i][j] = -lbm.u0 * std::cos(x) * std::sin(y);
      lbm.v[i][j] = lbm.u0 * std::sin(x) * std::cos(y);
    }
  }

  lbm.initializeDistributionFunction();

  std::cout << "finish initializing" << std::endl;

}


void simulateCavity2D(Parameters params, std::string dir, int idx, bool uq)
{

  lbm lbm(dir);
  lbm.UnitConverterFromResolutionAndRelaxationTime(params);

  setGeometry(lbm, params);

  initialize(lbm);

  std::cout << "start iteration" << std::endl;
  double td = 1.0 / (lbm.physViscosity * (lbm.dx * lbm.dx * 2.0));
  std::cout << "td: " << td << std::endl;
  int count = 0;
  //std::clock_t c_start = std::clock();
  double start = omp_get_wtime();
  double start_0 = start;
  //std::clock_t c_end = std::clock();
  double end = omp_get_wtime();
  double err = 0.;

  double t = 0.0, t0, t1;
  
  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;
  double error = 1.0;
  double uNorm = 1.0;

  #pragma omp parallel 
    while(error > 0.0001) {
      //for (int i = 1; i < 3; ++i) {
      lbm.collision();
      //sglbm.boundary();
      lbm.streaming();
      lbm.reconstruction(); 
      #pragma omp single
      {
        count++;
        if (count % 1000 == 0){
          error = calcError(lbm, uNorm);
          end = omp_get_wtime();
          std::cout << "iter: " << count << " " << "CPI time used: " << end - start << "s" << "\t" << "err: " << error << std::endl;
          start = end;
          t = 0.0;
          if ( count % 10000 == 0 && uq == false ) {
            lbm.output(dir, count, uq);
          }
        }
      }
  
    }

  std::cout << "total MCS time used: " <<  end - start_0 << std::endl;
  end = omp_get_wtime();
  count = int(td * 0.5) - 1;
  if (uq) {
    lbm.output(dir, idx, uq);
    velocity_All(dir, lbm, idx, uq);
  } else {
    lbm.output(dir, count, uq);
    velocity_central(dir, lbm, (lbm.nx+1)/2, (lbm.ny+1)/2);
    velocity_All(dir, lbm, 0, uq);
  }
    
}
