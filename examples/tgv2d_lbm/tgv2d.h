#include <cmath>
#include <random>
#include "../../src/lbm.h"
#include "../../src/postprocessing_lbm.h"

double calc_tke_error(lbm lbm, int count) {
  
  double tke = 0.0;
  double tkeAna = 0.0;
  totalKineticEnergy(lbm, tke, tkeAna, count);
  return std::abs((tke - tkeAna) / tkeAna);
}

void setGeometry(lbm& lbm, Parameters params) {
  std::cout << "start setting geometry" << std::endl;

  lbm.nx = lbm.N;
  lbm.ny = lbm.N;

  lbm.material = std::vector<std::vector<int>>(lbm.nx, std::vector<int>(lbm.ny, 1));

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

}


void simulateTGV2D(Parameters params, std::string dir, int idx, bool uq)
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

  #pragma omp parallel 
    for (int i = 1; i < int(td * 0.5); ++i) {
      //for (int i = 1; i < 3; ++i) {
      lbm.collision();
      //sglbm.boundary();
      lbm.streaming();
      lbm.reconstruction(); 
  
    }

  std::cout << "total MCS time used: " <<  end - start_0 << std::endl;
  end = omp_get_wtime();
  count = int(td * 0.5) - 1;
  if (uq) {
    lbm.output(dir, idx, uq);
    velocity_All(dir, lbm, idx, uq);
  } else {
    lbm.output(dir, count, uq);
    velocity_All(dir, lbm, 0, uq);
  }
  
  outputTKE(dir, lbm, count, end - start_0, uq);
  
}
