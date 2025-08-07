#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include "main_new.h"
#include "Schemes.h"
//instances of the "Fullerene" class have all the physical characteristics of a fullerene: 
//the speed of the center of mass, 
//coordinates, kinetic momentum, rotation speed
class Fullerene {
    int n1=2;
    public:
    std::vector <double> Xk, Yk, Zk;
    std::vector<double> tenz_ihertia;
    double m= 1.992646 * pow(10, -26);
    std::vector <double> xk, yk, zk;
    double uc, vc ,wc,yc, xc ,zc;
    double p, q, r;
    double Kx, Ky, Kz;
    Fullerene() {
      Xk.resize(n1), Yk.resize(n1), Zk.resize(n1);
      xk.resize(60), yk.resize(60), zk.resize(60);
        }
    Fullerene(std::fstream &f1,
              std::fstream &f2,
              std::fstream &f3,
              double u0,
              double v0,
              double w0,
              double p0,
              double q0,
              double r0) {
        Xk.resize(n1), Yk.resize(n1), Zk.resize(n1);
        xk.resize(60), yk.resize(60), zk.resize(60);
 for (int i=0;i<xk.size();i++) {
           f1>>xk.at(i);
           f2>>yk.at(i);
           f3>>zk.at(i);
           xk.at(i)*=0.1;
           yk.at(i)*=0.1;//перевод в нм
           zk.at(i)*=0.1;}
           uc=u0, vc=v0, wc=w0;
           p=p0, q=q0, r=r0;
    }
   Fullerene(std::fstream &f1,
             std::fstream &f2,
             std::fstream &f3,
             double d,
             double u0,
             double v0,
             double w0,
             double p0,
             double q0,
             double r0) {
        Xk.resize(n1), Yk.resize(n1), Zk.resize(n1);
        xk.resize(60), yk.resize(60), zk.resize(60);

        for (int i=0;i<xk.size();i++) {
           f1>>xk.at(i);
           f2>>yk.at(i);
           f3>>zk.at(i);

           xk.at(i)*=0.1;
           yk.at(i)*=0.1;//перевод в нм
           zk.at(i)*=0.1;

           xk.at(i)-=d;
           yk.at(i)-=d;
           zk.at(i)-=d;
        }
        uc=u0, vc=v0, wc=w0;
        p=p0, q=q0, r=r0;
    }

     void CalculateToTenz_ihertia(){
    for (int i=0;i<tenz_ihertia.size();i++)
    {
        tenz_ihertia.at(0)+=(pow(yk.at(i)-yc,2) + pow(zk.at(i)-zc,2));//A
        tenz_ihertia.at(1)+=(pow(xk.at(i)-xc,2) + pow(zk.at(i)-zc,2));//B
        tenz_ihertia.at(2)+=(pow(yk.at(i)-yc,2) + pow(xk.at(i)-xc,2));//C

        tenz_ihertia.at(3)-=(yk.at(i)-yc) * (zk.at(i)- zc);//D
        tenz_ihertia.at(4)-=(xk.at(i)-xc) * (zk.at(i)- zc);//E
        tenz_ihertia.at(5)-=(xk.at(i)-xc) * (yk.at(i)- yc);//F
    }
   tenz_ihertia.at(0)*=m; tenz_ihertia.at(1)*=m; tenz_ihertia.at(2)*=m; tenz_ihertia.at(3)*=m; 
   tenz_ihertia.at(4)*=m; tenz_ihertia.at(5)*=m; 
 }
void CalculateToKineticMoment() {
  Kx = tenz_ihertia.at(0) * p + tenz_ihertia.at(5) * q +tenz_ihertia.at(4) * r;
  Ky = tenz_ihertia.at(5) * p +tenz_ihertia.at(1) * q + tenz_ihertia.at(3) * r;
  Kz = tenz_ihertia.at(4) * p + tenz_ihertia.at(3) * q + tenz_ihertia.at(2) * r;
 }


};
//the class "Fullerene_new" is needed to obtain shifted coordinates 
//for calculating the coefficients of the R-K method
class Fullerene_new: public Fullerene {
 public:
 Fullerene_new(std::fstream &f1,
              std::fstream &f2,
              std::fstream &f3,
              double u0,
              double v0,
              double w0,
              double p0,
              double q0,
              double r0): Fullerene(f1,f2,f3,u0,v0,w0,p0,q0,r0) {}
Fullerene_new(std::fstream &f1,
              std::fstream &f2,
              std::fstream &f3,
              double d,
              double u0,
              double v0,
              double w0,
              double p0,
              double q0,
              double r0): Fullerene(f1,f2,f3,d,u0,v0,w0,p0,q0,r0) {}
//coordinate and velocity shifts are needed to c
//alculate the coefficients of the method (k1x, k2x,..., k1u ,..., k4u)
 double sdvig(Schemes &S,Fullerene &F1, int l){
    //Stage 1 of the method R-K
    uc = F1.uc + S.h * S.kuc.at(l);
    vc = F1.vc + S.h * S.kvc.at(l);
    wc = F1.wc + S.h * S.kwc.at(l);

    xc = F1.xc + S.h * S.kxc.at(l);
    yc = F1.yc + S.h * S.kyc.at(l);
    zc = F1.zc + S.h * S.kzc.at(l);


 }


};


