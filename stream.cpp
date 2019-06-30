//
// Created by Priyam Dalmia on 14-05-2019.
//
#include<iostream>
#include <vector>
#include <math.h>
using namespace std;
int nx = 128;
int ny = 128;

void stream(double fn[][130][9]) {

    for (int i = 0; i < ny+2; i++) {
        for (int j = nx; j >= 0; j--) {
            fn[i][j + 1][1] = fn[i][j][1];
        }
    }
///////////////////////////////////////////
    for (int i = 0; i < ny+2; i++) {
        for (int j = 1; j < nx+2; j++) {
            fn[i][j - 1][3] = fn[i][j][3];
        }
    }
///////////////////////////////////////////////////

    for (int j = 0; j < nx+2; j++) {
        for (int i = ny; i >= 0; i--) {
            fn[i + 1][j][4] = fn[i][j][4];
        }
    }
////////////////////////////////////////////////
    for (int j = 0; j < nx+2; j++) {
        for (int i = 1; i < ny+2; i++) {
            fn[i - 1][j][2] = fn[i][j][2];
        }
    }
/////////////////////////////////////////////////////
    for (int j = 0; j < nx+1; j++) {
        int k = j;
        for (int i = ny; i >= ny - j; i--) {
            fn[i + 1][k + 1][5] = fn[i][k][5];
            k--;
        }
    }

    for (int i = ny-1; i >= 0; i--) {
        int k = i;
        for (int j = nx; j >= nx - i; j--) {
            fn[k + 1][j + 1][5] = fn[k][j][5];
            k--;
        }
    }
///////////////////////////////////////////////////////
    for (int j = 0; j < nx+1; j++) {
        int k = j;
        for (int i = 1; i <= j + 1; i++) {
            fn[i - 1][k + 1][6] = fn[i][k][6];
            k--;
        }
    }
    for (int i = 2; i <= ny+1; i++) {
        int k = i;
        for (int j = nx; j >= i - 1; j--) {
            fn[k - 1][j + 1][6] = fn[k][j][6];
            k++;
        }
    }

//////////////////////////////////////////////////
    for (int j = 1; j <= nx+1; j++) {
        int k = j;
        for (int i = 1; i <= ny+1 - j + 1; j++) {
            fn[i - 1][k - 1][7] = fn[i][k][7];
            k++;
        }
    }
    for (int i = 2; i <= ny+1; i++) {
        int k = i;
        for (int j = 1; j <= nx+1 - i + 1; j++) {
            fn[k - 1][j - 1][7] = fn[k][j][7];
            k++;
        }
    }
//////////////////////////////////////////////////

    for (int j = 1; j <= nx+1; j++) {
        int k = j;
        for (int i = ny; i >= j - 1; i--) {
            fn[i + 1][k - 1][8] = fn[i][k][8];
            k++;
        }
    }
    for (int i = 0; i < ny; i++) {
        int k = i;
        for (int j = 1; j <= i + 1; j++) {
            fn[k + 1][j - 1][8] = fn[k][j][8];
            k--;
        }
    }

}

void applybc(double fn[][130][9],double f0[],double f1[]) {

    double sumLeft   = f0[3] + f0[7] + f0[8];
    double sumTop    = f1[4] + f1[5] + f1[8];
    double sumRight  = f0[1] + f0[5] + f0[6];
    double sumBottom = f0[7] + f0[6] + f0[2];
    
    int i = ny, l = 1;
    
    for (int j = 1; j < nx+1; j++) {
        
        double sum  = fn[i+1][j][4] + fn[i+1][j][5] + fn[i+1][j][8];
        fn[i][j][2] = (sum / sumTop) * f1[2];
        fn[i][j][6] = (sum / sumTop) * f1[6];
        fn[i][j][7] = (sum / sumTop) * f1[7];

        double suml = fn[l-1][j][2] + fn[l-1][j][6] + fn[l-1][j][7];
        fn[l][j][4] = (suml / sumBottom) * f0[4];
        fn[l][j][5] = (suml / sumBottom) * f0[5];
        fn[l][j][8] = (suml / sumBottom) * f0[8];
    }
    
    int j = 1;
    l = nx;
    for (i = 1; i < ny+1; i++) {
        double sum = fn[i][j-1][3] + fn[i][j-1][7] + fn[i][j-1][8];
        double suml = fn[i][l+1][1] + fn[i][l+1][5] + fn[i][l+1][6];
        fn[i][j][1] = (sum / sumLeft) * f0[1];
        fn[i][j][5] = (sum / sumLeft) * f0[5];
        fn[i][j][6] = (sum / sumLeft) * f0[6];
        fn[i][l][3] = (suml / sumRight) * f0[3];
        fn[i][l][7] = (suml / sumRight) * f0[7];
        fn[i][l][8] = (suml / sumRight) * f0[8];
    }
}

void getnew(double rho[][130],double u[][130],double v[][130],double fn[][130][9]) {

    for (int i = 1; i < ny+1; i++) {
        for (int j = 1; j < nx+1; j++) {
            rho[i][j] = 0;
            for (int k = 0; k < 9; k++) {
                rho[i][j] = rho[i][j] + fn[i][j][k];
            }
            u[i][j] = (fn[i][j][1] - fn[i][j][3] + fn[i][j][5] + fn[i][j][6] - fn[i][j][7] - fn[i][j][8])/rho[i][j];
            v[i][j] = (-fn[i][j][2] + fn[i][j][4] + fn[i][j][5] - fn[i][j][6] - fn[i][j][7] + fn[i][j][8])/rho[i][j];
            }
        }
}

void getFeq(double rho, double u, double v, double (&fEq)[9])
{
    double a,b,c;
    double f[9] ;
    //a,b,c are functions of the lagrange multipliers
    b = ((2.0 * u) + sqrt(1.0 + (3.0 * u * u))) / (1.0 - u);
    c = ((2.0 * v) + sqrt(1.0 + (3.0 * v * v))) / (1.0 - v);
    a = 36.0 / ((4.0 + b + (1 / b)) * (4.0 + c + (1 / c)));
    // cout<<a <<" "<<b<<" "<<c;
    // f-eq defined as the following
    f[0] = rho *(16.0 * a / 36.0);
    f[1] = rho * a * b * 4.0 / 36.0;
    f[2] = rho *(a / c) * 4.0 / 36.0;
    f[3] = rho *(a / b) * 4.0 / 36.0;
    f[4] = rho *a * c * 4.0 / 36.0;
    f[5] = rho *a * b * c / 36.0;
    f[6] = rho *a * b / (c * 36.0);
    f[7] = rho *a / (b * c * 36.0);
    f[8] = rho *a * c / (b * 36.0);
    // return &f[0];
    
//    double temp(0.0);
//    for(int dv=0;dv<9;dv++){
//
//      fEq[dv] = f[dv];
////     std::cout << fEq[dv] << "  ";
//    }
      
}

void getGlobalMass(double fn[][130][9] ) {
    double temp(0.0);
    for (int i = 1; i < ny+1; i++) {
        for (int j = 1; j < nx+1; j++) {
            for (int k = 0; k < 9; k++) {
                temp+= fn[i][j][k] ;
            }
        }
    }
    std::cout << "mass="<<temp<<"\n";
}

void collision(double b, double rho[][130],double u[][130],double v[][130],double fn[][130][9]) {
  
    double feqm[9],temp;
    for (int i = 1; i < ny+1; i++)
    {
      for (int j = 1; j < nx+1; j++)
      {
	  getFeq(rho[i][j], u[i][j], v[i][j], feqm);

	  for (int k = 0; k < 9; k++) {
	  
                fn[i][j][k] = fn[i][j][k] + 2.0*b*( feqm[k] - fn[i][j][k] );
            }

      }
    }
}

void bcwall(double fn[][130][9]) {
    for (int j = 0; j < nx+2; j++) {
        for (int k = 0; k < 9; k++) {
            fn[0][j][k] = fn[1][j][k];
            fn[129][j][k] = fn[128][j][k];
        }

    }
    for (int j = 0; j < ny+2; j++) {
        for (int k = 0; k < 9; k++) {
            fn[j][129][k] = fn[j][128][k];
            fn[j][0][k] = fn[j][1][k];
        }

    }
}


void bcperiodic(double fn[][130][9]) {
    for (int j = 0; j < nx+2; j++) {
        for (int k = 0; k < 9; k++) {
            fn[129][j][k] = fn[1][j][k];
            fn[0][j][k] = fn[128][j][k];
        }

    }
    for (int j = 0; j < ny+2; j++) {
        for (int k = 0; k < 9; k++) {
            fn[j][0][k] = fn[j][128][k];
            fn[j][129][k] = fn[j][1][k];
        }

    }
}