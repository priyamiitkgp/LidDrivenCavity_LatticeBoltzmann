#include<fstream>
#include<iomanip>
#include<time.h>
#include<string.h>
#include<stdio.h>

#include <iostream>
#include <vector>
#include<math.h>
#include<stdlib.h>
#include <cstdlib>
#include<cmath>

#include"stream.cpp"


using namespace std;

int main() {
    //Declarations

    double fn[130][130][9], rho[130][130], u[130][130], v[130][130];
    double f[9];
    double b, tau, ux;
    double Ma = 0.1;
    double Re = 100.0;
    
    tau = (Ma * sqrt(3) * 128) / Re;
    b = 1.0 / (1.0 + (2 * tau));
    ux = Ma / sqrt(3);
    
    int t = 50 * nx / ux;
    double f1[9], f0[9];
    getFeq(1, ux, 0, f1);
    getFeq(1, 0, 0, f0);
    
    std::cout << "beta="<<b<<"\n";

    //initialize f for first iteration
    for (int i = 1; i < ny+1; i++) {
        for (int j = 1; j < nx+1; j++) {
            float rp;
            rp = rand() % 1000;
            double r = rp / 1000;
            rho[i][j] = 1;
            u[i][j] = 0;
            v[i][j] = 0;
            getFeq(rho[i][j], u[i][j], v[i][j], f);
            for (int k = 0; k < 9; k++) {
                fn[i][j][k] = f[k];
            }
        }
    }

    std::cout << "t=0"<<"  ";
    getGlobalMass(fn);

    //Time loop
    for (int h = 0; h < t; h++) {
 
        getnew(rho,u,v,fn);
 	    collision(b,rho,u,v,fn);
	    bcperiodic(fn);
        bcwall(fn);
        stream(fn);
        applybc(fn,f0,f1);

	if(h%10==0){
	  std::cout << "t="<<h<<"  ";
	  getGlobalMass(fn);
	}
    }
//    for (int i=1;i<129;i++){
//        cout<<" "<<u[i][64];
//    }
//for (int j=1;j<129;j++){
  //  cout<<" "<<v[64][j];
//}

}