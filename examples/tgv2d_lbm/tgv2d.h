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


void initialize(lbm& lbm, std::vector<double> delta) {

  lbm.prepareLattice();

  for (int i = 0; i < lbm.nx; ++i) {
    for (int j = 0; j < lbm.ny; ++j) {
      double x = i * lbm.dx;
      double y = j * lbm.dx;

      double perturbation = 0.;

      perturbation += delta[0] * std::sin(2*x) * std::sin(2*y);
      perturbation += delta[1] * std::sin(2*x) * std::cos(2*y);
      perturbation += delta[2] * std::cos(2*x) * std::sin(2*y);
      perturbation += delta[3] * std::cos(2*x) * std::cos(2*y);

              
      lbm.rho[i][j] = 1.0 - 1.5 * lbm.u0 * lbm.u0 * std::cos(x + y) * std::cos(x - y);
      lbm.u[i][j] = -lbm.u0 * (1 + 0.25 * perturbation) * std::cos(x) * std::sin(y);
      lbm.v[i][j] =  lbm.u0 * (1 + 0.25 * perturbation) * std::sin(x) * std::cos(y);
    }
  }

  lbm.initializeDistributionFunction();

}


void simulateTGV2D(Parameters params, std::string dir, int idx, bool uq)
{

  lbm lbm(dir);
  lbm.UnitConverterFromResolutionAndRelaxationTime(params);

  setGeometry(lbm, params);

  std::vector<double> delta(params.polynomialType.size(), 0.0);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::cout << "delta: ";
  for (int i = 0; i < params.polynomialType.size(); ++i) {
    std::uniform_real_distribution<> dis(params.parameter1[i], params.parameter2[i]);
    delta[i] = dis(gen);
    std::cout << delta[i] << "\t";
  }
  std::cout << std::endl;

  initialize(lbm, delta);

  std::cout << "start iteration" << std::endl;
  double td = 1.0 / (lbm.physViscosity * (lbm.dx * lbm.dx * 2.0));
  std::cout << "td: " << td << std::endl;
  int count = 0;
  
  double start = omp_get_wtime();
  double start_0 = start;
  double err = 0.;
  
  size_t cores = omp_get_num_procs();
  omp_set_dynamic(0);
  omp_set_num_threads(cores);
  std::cout << "num Threads: " << cores << std::endl;

  #pragma omp parallel 
    for (int i = 1; i < int(td * 0.5); ++i) {
      lbm.collision();
      lbm.streaming();
      lbm.reconstruction();   
    }

  
  double end = omp_get_wtime();
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
