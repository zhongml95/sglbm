#include <iostream>
#include <cmath>
//#include "src/polynomial.h"
#include "src/sglbm.h"


int main( int argc, char* argv[] )
{
    unsigned int order = 1;
    unsigned int nq = 100;
    int resolution = 128;
    double parameter1 = 0.9;
    double parameter2 = 1.1;
    double L = 2.0 * M_PI;
    double lx = 2.0 * M_PI;
    double ly = 2.0 * M_PI;
    double dx = L / resolution;
    double dy = L / resolution;
    int nx = int(lx / dx);//+1;
    int ny = int(ly / dy);//+1;

    double tau = 0.5125663706143592;
    double Re = 15;
    double physVelocity = 0.01;
    double nu = physVelocity*2*M_PI/Re;
    std::vector<std::vector<int>> material(resolution+1, std::vector<int>(resolution+1, 1));
    
    std::string dir = "./data/tgv/t5/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(nx) + "/";
    std::string dirAna = "./data/tgv/t5/Nr" + std::to_string(order) + "Nq" + std::to_string(nq) + "N" + std::to_string(nx) + "/final/";

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
        sglbm.iteration();



    return 0;
}