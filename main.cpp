#include <iostream>
#include <cmath>
//#include "src/polynomial.h"
#include "src/sglbm.h"


int main( int argc, char* argv[] )
{
    unsigned int order = 4;
    unsigned int nq = 9;
    int resolution = 32;
    double parameter1 = 0.8;
    double parameter2 = 1.2;
    double L = 2.0 * M_PI;
    double lx = 2.0 * M_PI;
    double ly = 2.0 * M_PI;

    double tau = 0.5125663706143592;
    double Re = 15;
    double physVelocity = 0.01;
    double nu = physVelocity*2*M_PI/Re;
    std::vector<std::vector<int>> material(resolution+1, std::vector<int>(resolution+1, 1));
    
    sglbm sglbm(nq, order, parameter1,parameter2);
    sglbm.setGeometry(L,resolution,lx,ly,material);
    sglbm.setFluid(physVelocity,nu,tau);
    sglbm.initialize();
    sglbm.iteration();


    return 0;
}