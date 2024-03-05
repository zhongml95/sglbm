#include "../../src/sglbm.h"

double calcError(sglbm sglbm, double&uNorm0){
  double error = 0.0;
  double uNorm1 = 0.0;
  for (int i = 0; i < sglbm.nx; ++i){
    for (int j = 0; j < sglbm.ny; j++){
      std::vector<double> uSlice(sglbm.op.order+1, 0.0);
      std::vector<double> vSlice(sglbm.op.order+1, 0.0);
      for (int alpha = 0; alpha < sglbm.op.order+1; ++alpha){
        uSlice[alpha] = sglbm.u[i][j][alpha];
        vSlice[alpha] = sglbm.v[i][j][alpha];
      }
      double uMean = 0.0;
      double vMean = 0.0;
      uMean = sglbm.op.mean(uSlice);
      vMean = sglbm.op.mean(vSlice);
      uNorm1 += (uMean * uMean + vMean * vMean);
    }
  }

  uNorm1 = std::sqrt(uNorm1);
  error = std::fabs(uNorm1 - uNorm0) / uNorm0;
  uNorm0 = uNorm1;
  return error;
}


int main( int argc, char* argv[] )
{
    Parameters params;
    
    // Call readParameters to populate the params instance
    readParameters("./parameters.dat", params);

    double dx = params.L / params.resolution;
    double dy = params.L / params.resolution;
    int nx = int(params.lx / dx)+1;
    int ny = int(params.ly / dy)+1;

    double physViscosity = params.physVelocity * params.L / params.Re;
    double tau = 0.5384;

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

    
    sglbm sglbm(dir, "cavity2d", params);
    sglbm.setGeometry(params.L, params.resolution, params.lx, params.ly, material);
    sglbm.setFluid(params.physVelocity, physViscosity, tau);
    sglbm.initialize();
    //sglbm.iteration();

    std::cout << "start iteration" << std::endl;
    int count = 1;
    //std::clock_t c_start = std::clock();
    double start = omp_get_wtime();
    //std::clock_t c_end = std::clock();
    double end = omp_get_wtime();

    double t = 0.0, t0, t1;

    sglbm.output(dir, 0, 0);


    size_t cores = omp_get_num_procs();
    omp_set_dynamic(0);
    omp_set_num_threads(8);

  double error = 1.0;
  double uNorm = 1.0;


#pragma omp parallel 
    //for (int i = 1; i < 100000; ++i) {
      while(error > 0.0001) {
      sglbm.collision();  // parallel for
      sglbm.boundary();
      sglbm.streaming();
      //}
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
            sglbm.output(dir, count, 0);
          }
        }
      }
    }
    sglbm.output(dir, count, 0);

    return 0;
}

