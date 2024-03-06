#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

// 格子数
const int imax = 200;
// ステップ数
const int nmax = 200;
// cfl数
const double cfl = 0.1;

// const double dx = ;
// const double dt = ;
const double dt_dx = 0.05;

double x[imax+1];
double u[imax+1];
double uic[imax+1];
double unew[imax+1];
double uexact[imax+1];
double uL[imax+1];
double uR[imax+1];

double minmod(double d1, double d2);
void reconst(double u[], double uL[], double uR[]);

int main(){
    // 計算格子
    for (int i=0; i<=imax; i++){
        x[i] = double(i);
    }
    // 初期条件
    for (int i=0; i<=imax; i++){
        if (x[i] < imax/2){
            uic[i] = 2.0;
        }
        else{
            uic[i] = 0.0;
        }
    }
    // 理論解
    double travel = double(nmax)*cfl;
    for (int i=0; i<=imax; i++){
        if (x[i] < imax/2+travel){
            uexact[i] = 2.0;
        }
        else{
            uexact[i] = 0.0;
        }
    }
    
    for (int i=0; i<=imax; i++){
        u[i] = uic[i];
    }

    for (int n=1; n<=nmax; n++){
        // セルを再構築
        reconst(u, uL, uR);
        // リーマン問題を解く
        for (int i=1; i<=imax-1; i++){
            // 
            double fR = 0.5*(uL[i+1]*uL[i+1]+uR[i]*uR[i])\
            -0.5*std::abs(uL[i+1]+uR[i])*0.5*(uL[i+1]-uR[i]);
            double fL = 0.5*(uL[i]*uL[i]+uR[i-1]*uR[i-1])\
            -0.5*std::abs(uL[i]+uR[i-1])*0.5*(uL[i]-uR[i-1]);
            unew[i] = u[i]-dt_dx*(fR-fL);
        }

        // 境界条件
        unew[0] = unew[1];
        unew[imax] = unew[imax-1];

        // 更新
        for (int i=0; i<=imax; i++){
            u[i] = unew[i];
        }
        if (n==int(nmax/2)){
            std::ofstream u_trans;
            u_trans.open("u_muscl_minmod_trans.csv", std::ios::trunc);

            for (int i=0; i<=imax; i++){
				u_trans << x[i] << "," << u[i] << std::endl;
			}
            u_trans.close();
        }
    }

    std::ofstream u_muscl;
    u_muscl.open("u_muscl_minmod_final.csv", std::ios::trunc);

    for (int i=0; i<=imax; i++){
        u_muscl << x[i] << "," << u[i] << "," << uic[i] << "," << uexact[i] << std::endl;
    }
    u_muscl.close();
}

// minmod

double minmod(double d1, double d2){
    if (d1*d2<=0){
        return 0.0;
    }
    else if(std::abs(d1) < std::abs(d2)){
        return d1;
    }
    else{
        return d2;
    }
}

// セルの再構築用関数
void reconst(double u[], double uL[], double uR[]){
    for (int i=1; i<=imax; i++){
        // delta-, delta+
        double dul = u[i]-u[i-1];
        double dur = u[i+1]-u[i];

        // 制限付き勾配
        double k = -1.0; // この値を変更すれば，計算精度が変わる
        double b = (3.0-k)/(1.0-k);
        double dul_lim = minmod(dur, b*dul);
        double dur_lim = minmod(dul, b*dur);

        uL[i] = u[i]-0.25*((1+k)*dul_lim+(1-k)*dur_lim);
        uR[i] = u[i]+0.25*((1-k)*dul_lim+(1+k)*dur_lim);
    }
    // 境界条件;
    uL[0] = uL[1];
    uR[0] = uR[1];

    uL[imax] = uL[imax-1];
    uR[imax] = uR[imax-1];
}