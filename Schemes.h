#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include "Fulleren.h"
//The "Schemes" class stores attributes of the numerical method
// for solving the system of equations of motion.
// Such attributes include all coefficients of this method, the time step.
class Schemes {
   double h = 1e-6;
   public:
    double const Bolcman_constant = 1.38 * pow(10, -23);
    const  int n1=60, n2=60;
    std:: vector <std:: vector <double>> kx;
    std:: vector <std:: vector <double>> ky;
    std:: vector <std:: vector <double>> kz;

    std:: vector <double> kxc, kyc, kzc;
    std:: vector <double> kuc, kvc, kwc;
    std:: vector <double> kKx, kKy, kKz;

    Schemes(){
     for (int i = 0; i < 4; ++i) {
       kx.push_back(std::vector<double>(60, 0.0));
       ky.push_back(std::vector<double>(60, 0.0));
       kz.push_back(std::vector<double>(60, 0.0));
    }
       kxc.resize(4), kyc.resize(4), kzc.resize(4);
       kuc.resize(4), kvc.resize(4), kwc.resize(4);
       kKx.resize(4), kKy.resize(4), kKz.resize(4);
    }
//calculation of coefficients kx, ky, kz. It is carried out in 
//the same way regardless of the stage number of the R-K method
double CalculateCoefficient(Fullerene_new &F1, int s) {
 for (int l=0;l<n1;l++) {
    kx.at(l).at(s) = F1.uc + F1.q * (F1.zk[l] - F1.zc) - F1.r * (F1.yk[l] - F1.yc);
    ky.at(l).at(s) = F1.vc + F1.r * (F1.xk[l] - F1.xc) - F1.p * (F1.zk[l] - F1.zc);
    kz.at(l).at(s) = F1.wc + F1.p * (F1.yk[l] - F1.yc) - F1.q * (F1.xk[l] - F1.xc);
    }
}

    ~ Schemes(){}

     
};