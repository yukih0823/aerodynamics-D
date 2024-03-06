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

double x[imax+1];
double u[imax+1];
double uic[imax+1];
double unew[imax+1];
double uexact[imax+1];

int main(){
	// 計算格子設定
	for (int i=0; i<=imax; i++){
		x[i] = double(i);
	}
	// 初期条件設定
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
	// 数値解
	for (int i=0; i<=imax; i++){
		u[i] = uic[i];
	}

	for (int n=1; n<=nmax; n++){
		for (int i=1; i<=imax-1; i++){
			unew[i] = u[i]-cfl*(u[i]-u[i-1]);
		}
		for (int i=1; i<=imax-1; i++){
			u[i] = unew[i];
		}
		if (n==int(nmax/2)){
			std::ofstream u_trans;
			u_trans.open("u_upw1_trans.csv", std::ios::trunc);

			for (int i=0; i<=imax; i++){
				u_trans << x[i] << "," << u[i] << std::endl;
			}
			u_trans.close();
		}
	}

	std::ofstream u_dat;
	u_dat.open("u_upw1_final.csv", std::ios::trunc);

	for (int i=0; i<=imax; i++){
		u_dat << x[i] << "," << u[i] << "," << uic[i] << "," << uexact[i] << std::endl;
 	}
	u_dat.close();
}