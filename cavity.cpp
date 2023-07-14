#include <iostream>
#include <cmath>
//#include "src/polynomial.h"
#include "src/sglbm.h"


int main( int argc, char* argv[] )
{
    unsigned int order = 7;
    unsigned int nq = 100;
    int resolution = 128;
    double parameter1 = 0.9;
    double parameter2 = 1.1;
    double L = 1.0;//2.0 * M_PI;
    double lx = 1.0;//2.0 * M_PI;
    double ly = 1.0;//2.0 * M_PI;

    double tau = 0.5384;
    double Re = 1000;
    double physVelocity = 0.2;
    double nu = 0.001;
    std::vector<std::vector<int>> material(resolution+1, std::vector<int>(resolution+1, 1));

    double dx = L / resolution;
    double dy = L / resolution;
    int nx = int(lx / dx);//+1;
    int ny = int(ly / dy);//+1;

    std::string dir = "./data/cavity2d/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(nx) + "/";
    std::string dirAna = "./data/cavity2d/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(nx) + "/final/";

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


    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            if (j == ny-1)
            {
                material[i][j] = 3;
            }
            if ((i == 0) || (i == nx-1) || (j == 0) || (j == ny-1)){
                material[i][j] = 2;
            }
        }
    }
    
    sglbm sglbm(dir, "cavity2d", nq, order, parameter1,parameter2);
    sglbm.setGeometry(L,resolution,lx,ly,material);
    sglbm.setFluid(physVelocity,nu,tau);
    sglbm.initialize();
    sglbm.iteration();



    return 0;
}