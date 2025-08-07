#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include "Fulleren.h"
#include "Schemes.h"
//Lennrad Jones potential constants
const double e1 = 12.5 * 1.38 * pow(10, -23), sg1 = 0.335;
//distance calculation function
double distance(double x1, double y1, double z1,double x2, double y2, double z2){
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}
//method for calculating projections of the principal vector of external forces
void Power(Fullerene &F1,Fullerene &F2, Schemes S1, Schemes S2){
    S1.kuc.at(0)=0.0;
    S1.kvc.at(0)=0.0;
    S1.kwc.at(0)=0.0;

    S2.kuc.at(0)=0.0;
    S2.kvc.at(0)=0.0;
    S2.kwc.at(0)=0.0;
//taking into account pairwise interaction, calculating van der Waals forces
 for (int s=0;s<n1;s++) {
    for (int l=0;l<n2;l++) {
        R = distance(F1.xk[s], F1.yk[s], F1.zk[s],F2.xk[l], F2.yk[l], F2.zk[l]);
        dR_ds = Pt(R,e1, sg1);
        F1.Xk.at(s)  += dR_ds * (F1.xk[s]  - F2.xk[l] ) / R;
        F1.Yk.at(s)  += dR_ds * (F1.yk[s]   - F2.yk[l] ) / R;
        F1.Zk.at(s)  += dR_ds * (F1.zk[s]   - F2.zk[l] ) / R;

        F2.Xk.at(l)  -= dR_ds * (F1.xk[s]  - F2.xk[l] )  / R;
        F2.Yk.at(l)  -= dR_ds * (F1.yk[s]  - F2.yk[l] ) / R;
        F2.Zk.at(l)  -= dR_ds * (F1.zk[s]  - F2.zk[l] ) / R;
    }
    S1.kuc.at(0)+=F1.Xk.at(s);
    S1.kvc.at(0)+=F1.Yk.at(s);
    S1.kwc.at(0)+=F1.Zk.at(s);
}
for (int l=0;l<n2;l++) {
      S2.kuc.at(0)+=F2.Xk.at(l);
      S2.kvc.at(0)+=F2.Yk.at(l);
      S2.kwc.at(0)+=F2.Zk.at(l);
    }
//calculation of coefficients of the k1u, k1v, k1w method 
//for both fullerenes
    S1.kuc.at(0)/=(60*F1.m);
    S1.kvc.at(0)/=(60*F1.m);
    S1.kwc.at(0)/=(60*F1.m);

    S2.kuc.at(0)/=(60*F2.m);
    S2.kvc.at(0)/=(60*F2.m);
    S2.kwc.at(0)/=(60*F2.m);
}

int main() {
int steps= pow(10,6);
std::ifstream f1("x_motion.txt");
std::ifstream f2("y_motion.txt");
std::ifstream f3("z_motion.txt");

//create the first fullerene by calling the first constructor
Fullerene F1(f1,f2,f3,100,150,-10,15,25,50);
Fullerene_new F1_1(f1,f2,f3,100,150,-10,15,25,50);

F1.xc=0.0, F1.yc=0.0, F1.zc=0.0;
F1_1.xc=0.0,F1_1.yc=0.0, F1_1.zc=0.0;
F1.CalculateToTenz_ihertia();
F1.CalculateToKineticMoment();

//create the second fullerene by calling the second constructor
Fullerene F2(f1,f2,f3,2.5,100,150,0,-15,25,50);
Fullerene_new F2_2(f1,f2,f3,2.5,100,150,0,-15,25,50);
F2.xc=-2.5, F2.yc=-2.5, F2.zc=-2.5;
F2_2.xc=-2.5,F2_2.yc=-2.5, F2_2.zc=-2.5;
F2.CalculateToTenz_ihertia();
F2.CalculateToKineticMoment();

Schemes RungeKutta4_1();// атрибуты разностной схемы для 1-го фуллерена 
Schemes RungeKutta4_2();

for (int i=0; i<steps; ++i){
// step R-K #1
RungeKutta4_1.CalculateCoefficient(F1,0);
RungeKutta4_2.CalculateCoefficient(F2,0);
//null
std::fill(F1.Xk.begin(), F1.Xk.end(), 0.0);
std::fill(F1.Yk.begin(), F1.Yk.end(), 0.0);
std::fill(F1.Zk.begin(), F1.Zk.end(), 0.0);

std::fill(F2.Xk.begin(), F2.Xk.end(), 0.0);
std::fill(F2.Yk.begin(), F2.Yk.end(), 0.0);
std::fill(F2.Zk.begin(), F2.Zk.end(), 0.0);

Power(F1, F2, RungeKutta4_1, RungeKutta4_2);
RungeKutta4_1.kxc.at(0) = F1.uc,
RungeKutta4_1.kyc.at(0) = F1.vc,
RungeKutta4_1.kzc.at(0) = F1.wc,

RungeKutta4_2.kxc.at(0) = F2.uc,
RungeKutta4_2.kyc.at(0) = F2.vc,
RungeKutta4_2.kzc.at(0) = F2.wc;

F1_1.sdvig(RungeKutta4_1, F1,0);
F2_2.sdvig(RungeKutta4_2, F2,0);

}


return 0;
}