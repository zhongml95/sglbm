#ifndef SGLBM_H
#define SGLBM_H
#include "polynomial.h"
#include <iostream>
#include <fstream> 

class sglbm : public LegendrePoly {
public:
    int N;
    double L;
    double dx;
    double dy;
    double lx;
    double ly;
    int nx;
    int ny;
    double physVelocity;
    double Re;
    double physViscosity;
    double tau;
    double dt;
    double conversionViscosity;
    double conversionVelocity;
    double conversionDensity;
    double conversionMass;
    double conversionForce;
    double u0;
    double omega0;

    double cs2 = 1.0/3.0;
    std::vector<int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    std::vector<int> cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    std::vector<double> w = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    std::vector<int> kinv = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    std::vector<std::vector<int>> material;
    std::vector<std::vector<std::vector<double>>> bouzidiQ;
    std::vector<std::vector<std::vector<double>>> u;
    std::vector<std::vector<std::vector<double>>> v;
    std::vector<std::vector<std::vector<double>>> rho;
    std::vector<std::vector<std::vector<double>>> omega;
    std::vector<std::vector<std::vector<std::vector<double>>>> f;
    std::vector<std::vector<std::vector<std::vector<double>>>> feq;

    sglbm(int _nq, int _n, double _parameter1, double _parameter2):LegendrePoly(_nq, _n,_parameter1,_parameter2)
    {

    }

    std::vector<double> find_intersection(std::vector<double> center, double radius, std::vector<int> start_point, std::vector<int> end_point)
    {
        double x0 = center[0];
        double y0 = center[1];

        double x1 = start_point[0];
        double y1 = start_point[1];

        double x2 = end_point[0];
        double y2 = end_point[1];

        std::vector<double> inp(2,0.0);
        if (radius == 0){
            return {x1, y1};
        }

        if (x1 == x2){
            inp[0] = x1;
            if (std::abs(radius) >= std::abs(x1 - x0)){
                double p1 = y0 - std::sqrt(radius*radius - (x1-x0)*(x1-x0));
                double p2 = y0 + std::sqrt(radius*radius - (x1-x0)*(x1-x0));
                if (std::max(y1,y2)>=p2){
                    inp[2] = p2;
                }
                if (std::min(y1,y2)<=p1){
                    inp[2] = p1;
                }
            }
        }
        else{
            double k = (y1 - y2) / (x1 - x2);
            double b0 = y1 - k * x1;
            
            double a = k*k + 1;
            double b = 2.0*k*(b0-y0) - 2.0*x0;
            double c = (b0 - y0)*(b0 - y0) + x0*x0 - radius*radius;
            double delta = b*b - 4*a*c;
            if (delta >= 0){                
                double p1x = (-1.0*b - std::sqrt(delta)) / (2*a);
                double p2x = (-1.0*b + std::sqrt(delta)) / (2*a);
                double p1y = k * p1x + b0;
                double p2y = k * p2x + b0;
                if (p1x >= std::min(x1,x2) && p1x <= std::max(x1,x2)){
                    inp[1] = p1x;
                    inp[2] = p1y;
                }
                else{
                    inp[1] = p2x;
                    inp[2] = p2y;
                }
            }
        }
        return inp;
    }

    void setGeometry(double _L, double _N, double _lx, double _ly, const std::vector<std::vector<int>>& _material){
        std::cout << "start setting geometry" << std::endl;
        L = _L;
        N = _N;
        dx = L/N;
        dy = L/N;
        lx = _lx;
        ly = _ly;
        nx = int(lx/dx);//+1;
        ny = int(ly/dy);//+1;
        material.resize(nx);
        for (int i = 0; i < nx; i++){
            material[i].resize(ny);
            for (int j = 0; j < ny; j++){
                material[i][j] = _material[i][j];
            }
        }
        std::cout << "finish setting geometry" << std::endl;

    }
    
    void setFluid(double _physVelocity, double _nu, double _tau)
    {
        std::cout << "set simulation parameter" << std::endl;
        physVelocity = _physVelocity;
        //Re = _Re;
        physViscosity = _nu;//physVelocity*L/Re;
        Re = physVelocity*L/physViscosity;
        tau = _tau;
        dt = (tau - 0.5) / 3.0 * (dx*dx) / physViscosity;
        conversionViscosity = dx*dx/dt;
        conversionVelocity = dx/dt;
        conversionDensity = 1.0;
        conversionMass = conversionDensity*dx*dx*dx;
        conversionForce = conversionMass * dx / dt/dt;
        u0 = physVelocity / conversionVelocity;
        omega0 = 1.0 / (3 * physViscosity/conversionViscosity + 0.5);
    }

    void setCircle(double centerX, double centerY, double radius)
    {
        bouzidiQ.resize(nx);
        for (int i = 0; i < nx; i++){
            bouzidiQ[i].resize(ny);
            for (int j = 0; j < ny; j++){
                bouzidiQ[i][j].resize(9);
            }
        }

        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                if (material[i][j] == 3){
                    for (int k = 1; k < 9; k++){
                        std::vector<double> intersection(2,0.0);
                        intersection = find_intersection({centerX/dx+1,centerY/dx+1},L/dx,{i,j},{i+cx[k],j+cy[k]});
                        bouzidiQ[i][j][k] = 1.0 - std::sqrt((i-intersection[1])*(i-intersection[1]) + (j-intersection[2])*(j-intersection[2])) / sqrt(cx[k]*cx[k]+cy[k]*cy[k]);
                    }
                }
            }
        }
    }

    void initialize()
    {
        std::cout << "starting initialization" << std::endl;

        std::cout << "resolution " << ny << std::endl;
        std::cout << "nx = " << nx << "\t" << "ny = " << ny << std::endl;
        std::cout << "tau = " << tau << std::endl;
        std::cout << "nu = " << physViscosity << std::endl;
        std::cout << "Re = " << Re << std::endl;
        std::cout << "u0 = " << u0 << std::endl;
        std::cout << "converstionVelocity = " << conversionVelocity << std::endl;
        std::cout << "converstionViscosity = " << conversionViscosity << std::endl;


        u.resize(nx);
        v.resize(nx);
        rho.resize(nx);
        omega.resize(nx);
        f.resize(nx);
        feq.resize(nx);
        for (int i = 0; i < nx; i++){            
            u[i].resize(ny);
            v[i].resize(ny);
            rho[i].resize(ny);
            omega[i].resize(ny);
            f[i].resize(ny);
            feq[i].resize(ny);
            for (int j = 0; j < ny; j++){
                u[i][j].resize(order+1);
                v[i][j].resize(order+1);
                rho[i][j].resize(order+1);
                omega[i][j].resize(order+1);
                for (int alpha = 0; alpha < order+1; alpha++){
                    u[i][j][alpha] = 0.0;
                    v[i][j][alpha] = 0.0;
                    rho[i][j][alpha] = 1.0;
                    omega[i][j][alpha] = omega0;
                }
                
                f[i][j].resize(9);
                feq[i][j].resize(9);
                for (int k = 0; k < 9; k++){
                    f[i][j][k].resize(order+1);
                    feq[i][j][k].resize(order+1);
                    for (int alpha = 0; alpha < order+1; alpha++){
                        f[i][j][k][alpha] = 0.0;
                        feq[i][j][k][alpha] = 0.0;
                    }
                }
            }
        }
        

        std::vector<double> uChaos(order+1,0.0);
        std::vector<double> vChaos(order+1,0.0);
        std::vector<double> rChaos(order+1,0.0);
        std::vector<double> omegaChaos(order+1,0.0);
        std::vector<double> chaos(2,0.0);
        rChaos[0] = 1.0;
        std::cout << "omega: " <<  1.0 / (3 * physViscosity*parameter2/conversionViscosity + 0.5) << std::endl;
        chaos = convert2affinePCE(1.0 / (3 * physViscosity*parameter2/conversionViscosity + 0.5), 1.0 / (3 * physViscosity*parameter1/conversionViscosity + 0.5));
        omegaChaos[0] = chaos[0];
        omegaChaos[1] = chaos[1];

        std::cout << "chaos: " << chaos[0] << "\t" << chaos[1] << std::endl;

        for (int i = 0; i < nx; i++){           
            for (int j = 0; j < ny; j++){    
                //TGV            
                double x = (i-1)*dx;
                double y = (j-1)*dy;
                rChaos[0] = 1.0 - 1.5*u0*u0*std::cos(x+y)*std::cos(x-y);
                uChaos[0] = -u0 * std::cos(x) * std::sin(y);
                vChaos[0] =  u0 * std::sin(x) * std::cos(y);

                for (int alpha = 0; alpha < order+1; alpha++){
                    u[i][j][alpha] = uChaos[alpha];
                    v[i][j][alpha] = vChaos[alpha];
                    rho[i][j][alpha] = rChaos[alpha];
                    omega[i][j][alpha] = omegaChaos[alpha];
                }
            }
        }

        for (int i = 0; i < nx; i++){           
            for (int j = 0; j < ny; j++){
                std::vector<double> rRan(nq,0.0);
                std::vector<double> uRan(nq,0.0);
                std::vector<double> vRan(nq,0.0);
                    
                std::vector<double> rSlice(order+1,0.0);
                std::vector<double> uSlice(order+1,0.0);
                std::vector<double> vSlice(order+1,0.0);
                for (int alpha = 0; alpha < order+1; alpha++){
                    rSlice[alpha] = rho[i][j][alpha];
                    uSlice[alpha] =   u[i][j][alpha];
                    vSlice[alpha] =   v[i][j][alpha];
                }

                rRan = evaluatePCE(rSlice);
                uRan = evaluatePCE(uSlice);
                vRan = evaluatePCE(vSlice);

                for (int k = 0; k < 9; k++){                    
                    std::vector<double> feqRan(nq,0.0);
                    for (int sample = 0; sample < nq; sample++){ 
                        feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
                    }

                    std::vector<double> feqSlice = ran2chaos(feqRan);
                    for (int alpha = 0; alpha < order+1; alpha++){
                        feq[i][j][k][alpha] = feqSlice[alpha];
                        //std::cout << feqSlice[alpha] << std::endl;
                        f[i][j][k][alpha] = feqSlice[alpha];
                    }
                }
            }
        }
    }

    void setTGV()
    {
        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                double x = (i-1)*dx;
                double y = (j-1)*dy;
                double rTGV = 1.0 - 1.5*u0*u0*std::cos(x+y)*std::cos(x-y);
                double uTGV = -u0 * std::cos(x) * std::sin(y);
                double vTGV = u0 * std::sin(x) * std::cos(y);
            }
        }
    }

    void collision()
    {
        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                for (int k = 0; k < 9; k++){
                    std::vector<double> omegaSlice(order+1,0.0);
                    std::vector<double> fSlice(order+1,0.0);
                    std::vector<double> feqSlice(order+1,0.0);
                    std::vector<double> QSlice(order+1,0.0);
                    for (int alpha = 0; alpha < order+1; alpha++){
                        omegaSlice[alpha] = omega[i][j][alpha];
                        fSlice[alpha] = f[i][j][k][alpha];
                        feqSlice[alpha] = feq[i][j][k][alpha];
                        //std::cout << feq[i][j][k][alpha] << std::endl;
                    }
                    QSlice = collisionTerm(fSlice, feqSlice, omegaSlice);

                    for (int alpha = 0; alpha < order+1; alpha++){
                        f[i][j][k][alpha] = f[i][j][k][alpha] + QSlice[alpha];
                        //std::cout << QSlice[alpha] << std::endl;
                    }
                }
            }
        }
    }

    std::vector<double> collisionTerm(std::vector<double> _f, std::vector<double> _feq, std::vector<double> _omega)
    {
        std::vector<double> Q(order+1,0.0);

        for (int i = 0; i < order+1; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < order+1; j++) {
                for (int k = 0; k < order+1; k++) {
                    sum1 += _omega[j] * _feq[k] * t3Product[j][k][i];
                    //std::cout << t3Product[j][k][i] << std::endl;
                    sum2 += _omega[j] * _f[k] * t3Product[j][k][i];
                }
            }

            Q[i] = (sum1 - sum2) / t2Product[i];
        }
        return Q;
    }

    double equilibrium(double _r, double _u, double _v, int _cx, int _cy, double _w)
    {
        double t1 = _u * _u + _v * _v;
        double t2 = _u * _cx + _v * _cy;
        return _r * _w * (1.0 + 3.0*t2 + 4.5*t2*t2 - 1.5*t1);
    }

    void streaming()
    {
        std::vector<std::vector<std::vector<std::vector<double>>>> ftmp(nx, std::vector<std::vector<std::vector<double>>>(ny, std::vector<std::vector<double>>(9, std::vector<double>(order + 1, 0.0))));
        int new_i = 0;
        int new_j = 0;
        
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < 9; k++) {

                    new_i = i + cx[k];
                    new_j = j + cy[k];

                    if (new_i == nx)
                        new_i = 0;

                    if (new_i == -1)
                        new_i = nx-1;

                    if (new_j == ny)
                        new_j = 0;

                    if (new_j == -1)
                        new_j = ny-1;

                    for (int alpha = 0; alpha < order + 1; alpha++) {
                        ftmp[new_i][new_j][k][alpha] = f[i][j][k][alpha];
                    }
                    
                }
            }
        }
        f.swap(ftmp);
        // Copy values back to f
        /*for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < 9; k++) {
                    for (int alpha = 0; alpha < order + 1; alpha++) {
                        f[i][j][k][alpha] = ftmp[i][j][k][alpha];
                    }
                }
            }
        }*/
    }

    void boundary()
    {

    }

    void reconstruction()
    {
        for (int i = 0; i < nx; i++){
            for (int j = 0; j< ny; j++){
                std::vector<double> rRan(nq,0.0);
                std::vector<double> uRan(nq,0.0);
                std::vector<double> vRan(nq,0.0);
                std::vector<double> ruRan(nq,0.0);
                std::vector<double> rvRan(nq,0.0);
                
                std::vector<double> rSlice(order+1,0.0);
                std::vector<double> uSlice(order+1,0.0);
                std::vector<double> vSlice(order+1,0.0);
                

                for (int alpha = 0; alpha < order+1; alpha++){
                    for (int k = 0; k < 9; k++){
                        rSlice[alpha] += f[i][j][k][alpha];
                    }
                    rho[i][j][alpha] = rSlice[alpha];
                }
                
                
                rRan = evaluatePCE(rSlice);

                if (material[i][j] == 1){
                    //std::cout << "check" << std::endl;
                    std::vector<double> ruSlice(order+1,0.0);
                    std::vector<double> rvSlice(order+1,0.0);
                    for (int alpha = 0; alpha < order+1; alpha++){
                        for (int k = 0; k < 9; k++){
                            ruSlice[alpha] += f[i][j][k][alpha] * cx[k];
                            rvSlice[alpha] += f[i][j][k][alpha] * cy[k];
                        }
                    }
                    //std::cout << "check" << std::endl;
                    ruRan = evaluatePCE(ruSlice);
                    rvRan = evaluatePCE(rvSlice);

                    for (int sample = 0; sample < nq; sample++){
                        uRan[sample] = ruRan[sample] / rRan[sample];
                        vRan[sample] = rvRan[sample] / rRan[sample];
                    }

                    uSlice = ran2chaos(uRan);
                    vSlice = ran2chaos(vRan);

                    for (int alpha = 0; alpha < order+1; alpha++){
                        u[i][j][alpha] = uSlice[alpha];
                        v[i][j][alpha] = vSlice[alpha];
                    }
                }

                for (int k = 0; k < 9; k++){                    
                    std::vector<double> feqRan(nq,0.0);
                    for (int sample = 0; sample < nq; sample++){ 
                        feqRan[sample] = equilibrium(rRan[sample], uRan[sample], vRan[sample], cx[k], cy[k], w[k]);
                    }

                    std::vector<double> feqSlice = ran2chaos(feqRan);
                    for (int alpha = 0; alpha < order+1; alpha++){
                        feq[i][j][k][alpha] = feqSlice[alpha];
                    }
                }
                
            }
        }
    }

    void output(int iter)
    {
        std::string filename = "./data/tgv_" + std::to_string(iter) + ".dat";
        std::ofstream outputFile(filename);
        if (!outputFile) {
            std::cerr << "Error opening the file: " << filename << std::endl;
            return;
        }

        outputFile << "VARIABLE = \"X\", \"Y\", \"U_MEAN\", \"U_STD\", \"V_MEAN\", \"V_STD\", \"P_MEAN\", \"P_STD\"\n";
        outputFile << "ZONE I = " << nx << ", J = " << ny << ", F = POINT\n";
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = static_cast<double>(i) / (nx - 1);
                double y = static_cast<double>(j) / (ny - 1);
                std::vector<double> rSlice(order+1,0.0);
                std::vector<double> uSlice(order+1,0.0);
                std::vector<double> vSlice(order+1,0.0);
                for (int alpha = 0; alpha < order+1; alpha++){
                    rSlice[alpha] = rho[i][j][alpha];
                    uSlice[alpha] = u[i][j][alpha];
                    vSlice[alpha] = v[i][j][alpha];
                }
                outputFile << x << "\t" << y << "\t" << mean(rSlice) << "\t" << mean(uSlice)*conversionVelocity << "\t" << mean(vSlice)*conversionVelocity << "\t" << std(rSlice) << "\t" << std(uSlice)*conversionVelocity << "\t" << std(vSlice)*conversionVelocity << "\n";
            }
        }
    }



    void iteration()
    {
        std::cout << "start iteration" << std::endl;
        double td = 1.0 / (physViscosity * (dx*dx + dy*dy));
        std::cout << "td: " << td << std::endl;
        for (int i = 0; i < int(td*0.5); i++){
            collision();            
            boundary();
            streaming();
            reconstruction();
            if (i%10 == 0){
                std::cout << "iter: " << i << std::endl;
                output(i);
            }
        }
    }
};

#endif // LBM_H