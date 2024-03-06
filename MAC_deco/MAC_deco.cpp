#include <iostream>
#include <fstream>
#include <cmath>
#include <vector> 
#include <string> 
#include <chrono>

// 凸型になるように物体を設定(もう少し計算格子を細かくしたほうがいいかも)

// Reynolds number
const double re = 70.0; 

// CFL number
const double cfl = 0.2; 

// SOR params
const double omegap = 1.00;
const int maxitp = 100;
const double errorp = 1e-4;

// No. of Time Steps
const int nlast = 5000;

// set x-grid params
const int mx = 401;
const int i_1 = 96; 
const int i_2 = 106;
const int i_3 = 116;
const double dx = 1.0/(i_2-i_1);

// set y-grid params
const int my = 201; 
const int j_1 = 86;
const int j_2 = 96;
const int j_3 = 106; 
const int j_4 = 116;
const double dy = 1.0/(j_2-j_1);

// クーラン条件
const double dt = cfl*std::fmin(dx, dy);
double simu_time = 0.0; // シミュレーション内での経過時間

using vd1 = std::vector<double>; // 1d matrix
using vd2 = std::vector<std::vector<double> >; // 2d matrix

vd1 x(mx+1), y(my+1);
vd2 u(mx+2, vd1(my+2)), v(mx+2, vd1(my+2)), p(mx+2, vd1(my+2));
vd2 rhs(mx+2, vd1(my+2));
vd2 omega(mx+2, vd1(my+2));

// サブ関数の前方宣言
void set_grid(vd1 & x, vd1 & y);
void intcnd(vd2 & u, vd2 & v, vd2 & p);
void bcforp(vd2 & p);
void bcforv(vd2 & u, vd2 & v);
void poiseq(vd2 & u, vd2 & v, vd2 & p);
void veloeq(vd2 & u, vd2 & v, vd2 & p);
void solve_flow(vd1 & x, vd1 & y, vd2 & u, vd2 & v, vd2 & p);

// メイン関数(処理を実行)
int main(){
    set_grid(x, y);
    intcnd(u, v, p);
    bcforp(p);
    bcforv(u, v);
    solve_flow(x, y, u, v, p);
    return 0;
}

/* サブ関数の定義 */

// 計算格子の設定
void set_grid(vd1 & x, vd1 & y){
    double icent = (i_1+i_2)/2;
    double jcent = (j_2+j_3)/2;
    for (int i=1; i<=mx; i++){
        x[i] = dx*double(i-icent);
    }
    for (int j=1; j<=my; j++){
        y[j] = dy*double(j-jcent);
    }
}

// 初期条件の設定
void intcnd(vd2 & u, vd2 & v, vd2 & p){
    for (int i=1; i<=mx; i++){
        for (int j=1; j<=my; j++){
            u[i][j] = 1.0; // t = 0.0で計算格子内に，物体が突如現れたと考える.
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            omega[i][j] = 0.0;
        }
    }
}

// 圧力境界条件の設定
void bcforp(vd2 & p){
    for (int j=1; j<=my; j++){
        p[1][j] = 0.0; // inflow condition(i=1)
        p[mx][j] = 0.0; // downstream condition(i=mx)
    }
    for (int i=1; i<=mx; i++){
        p[i][1] = 0.0; // bottom condition(j=1)
        p[i][my] = 0.0; // bottom condition(j=my)
    }
    // wall condition
    // まずは隅を決定
    p[i_2][j_1] = p[i_2-1][j_1-1];
    p[i_2][j_2] = p[i_2-1][j_2-1];
    p[i_1][j_2] = p[i_1-1][j_2-1];
    p[i_1][j_3] = p[i_1-1][j_3+1];
    p[i_2][j_3] = p[i_2-1][j_3+1];
    p[i_2][j_4] = p[i_2-1][j_4+1];
    p[i_3][j_4] = p[i_3+1][j_4+1];
    p[i_3][j_1] = p[i_3+1][j_1-1];
    // 次に辺を決定
    // 横
    for (int j=j_2+1; j<=j_3-1; j++){
        p[i_1][j] = p[i_1-1][j];
    }
    for (int j=j_1+1; j<=j_2-1; j++){
        p[i_2][j] = p[i_2-1][j];
    }
    for (int j=j_3+1; j<=j_4-1; j++){
        p[i_2][j] = p[i_2-1][j];
    }
    for (int j=j_1+1; j<=j_4-1; j++){
        p[i_3][j] = p[i_3+1][j];
    }
    // 縦
    for (int i=i_1+1; i<=i_2-1; i++){
        p[i][j_2] = p[i][j_2-1];
        p[i][j_3] = p[i][j_3+1];
    }
    for (int i=i_2+1; i<=i_3-1; i++){
        p[i][j_1] = p[i][j_1-1];
        p[i][j_4] = p[i][j_4+1];
    }
}

// 速度境界条件の設定
void bcforv(vd2 & u, vd2 &v){
    for (int j=1; j<=my; j++){
        // inflow condition(i=1)
        u[1][j] = 1.0;
        v[1][j] = 0.0;
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        // downstream condition(i=mx)
        u[mx][j] = 2.0*u[mx-1][j]-u[mx-2][j];
        v[mx][j] = 2.0*v[mx-1][j]-v[mx-2][j];
        u[mx+1][j] = 2.0*u[mx][j]-u[mx-1][j];
        v[mx+1][j] = 2.0*v[mx][j]-v[mx-1][j];
    }
    for (int i=1; i<=mx; i++){
        // bottom condition(j=1)
        u[i][1] = 2.0*u[i][2]-u[i][3];
        v[i][1] = 2.0*v[i][2]-v[i][3];
        u[i][0] = 2.0*u[i][1]-u[i][2];
        v[i][0] = 2.0*v[i][1]-v[i][2];
        // bottom condition(j=my)
        u[i][my] = 2.0*u[i][my-1]-u[i][my-2];
        v[i][my] = 2.0*v[i][my-1]-v[i][my-2];
        u[i][my+1] = 2.0*u[i][my]-u[i][my-1];
        v[i][my+1] = 2.0*v[i][my]-v[i][my-1];
    }
    // wall condition
    for (int i=i_1; i<=i_2; i++){
        for (int j=j_2; j<=j_3; j++){
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
    for (int i=i_2; i<=i_3; i++){
        for (int j=j_1; j<=j_4; j++){
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

// 圧力場を解く(緩和法によるポアソン方程式の解)
void poiseq(vd2 & u, vd2 & v, vd2 & p){
    // vd2 rhs(mx+2, vd1(my+2)); // rhsをreturnしてステップごとの変化をみたい
    // ポアソン方程式の右辺をまず計算
    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            double ux = (u[i+1][j]-u[i-1][j])/(2.0*dx);
            double uy = (u[i][j+1]-u[i][j-1])/(2.0*dy);
            double vx = (v[i+1][j]-v[i-1][j])/(2.0*dx);
            double vy = (v[i][j+1]-v[i][j-1])/(2.0*dy);
            rhs[i][j] = (ux+vy)/dt-(ux*ux+2.0*uy*vx+vy*vy);
        }
    }

    // iterations
    for (int itr=1; itr<=maxitp; itr++){
        double res = 0.0;
        for (int i=2; i<=mx-1; i++){
            for (int j=2; j<=my-1; j++){
                // 物体内部は計算しない
                if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                    continue;
                }
                else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                    continue;
                }
                double dp = (p[i+1][j]+p[i-1][j])/(dx*dx)+(p[i][j+1]+p[i][j-1])/(dy*dy)-rhs[i][j];
                dp = dp/(2.0/(dx*dx)+2.0/(dy*dy))-p[i][j]; // p[i][j]の修正すべき量
                res = res+dp*dp;
                p[i][j] = p[i][j]+omegap*dp; // 修正すべき量に係数をかけて修正
            }
        }
        bcforp(p);
        res = std::sqrt(res/double(mx*my));
        if (res < errorp){
            break;
        }
    }
}

// 速度場を解く(Kawamuraスキームの適用)
void veloeq(vd2 & u, vd2 & v, vd2 & p){
    vd2 urhs(mx+2, vd1(my+2)), vrhs(mx+2, vd1(my+2)); // 圧力勾配, 粘性項，移流項などの速度変化の要因を計算
    // 圧力勾配
    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            urhs[i][j] = -(p[i+1][j]-p[i-1][j])/(2.0*dx);
            vrhs[i][j] = -(p[i][j+1]-p[i][j-1])/(2.0*dy);
        }
    }
    // 粘性項
    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            urhs[i][j] =
                        urhs[i][j]+
                        (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/(re*dx*dx)+
                        (u[i][j+1]-2.0*u[i][j]+u[i][j-1])/(re*dy*dy);
            vrhs[i][j] =
                        vrhs[i][j]+
                        (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/(re*dx*dx)+
                        (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/(re*dy*dy);
        }
    }
    // advection term in x-direction(kawamuraスキームの計算に必要)
    for (int j=j_2+1; j<=j_3-1; j++){
        u[i_1+1][j] = 2.0*u[i_1][j]-u[i_1-1][j];
        v[i_1+1][j] = 2.0*v[i_1][j]-v[i_1-1][j];
    }
    for (int j=j_1+1; j<=j_2-1; j++){
        u[i_2+1][j] = 2.0*u[i_2][j]-u[i_2-1][j];
        v[i_2+1][j] = 2.0*v[i_2][j]-v[i_2-1][j];
    }
    for (int j=j_3+1; j<=j_4-1; j++){
        u[i_2+1][j] = 2.0*u[i_2][j]-u[i_2-1][j];
        v[i_2+1][j] = 2.0*v[i_2][j]-v[i_2-1][j];
    }
    for (int j=j_1+1; j<=j_4-1; j++){
        u[i_3-1][j] = 2.0*u[i_3][j]-u[i_3+1][j];
        v[i_3-1][j] = 2.0*v[i_3][j]-v[i_3+1][j];
    }

    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            urhs[i][j] = urhs[i][j]-u[i][j]*(-u[i+2][j]+8.0*(u[i+1][j]-u[i-1][j])+u[i-2][j])/(12.0*dx)
                         -abs(u[i][j])*(u[i+2][j]-4.0*u[i+1][j]+6.0*u[i][j]-4.0*u[i-1][j]+u[i-2][j])/(4.0*dx);
			vrhs[i][j] = vrhs[i][j]-u[i][j]*(-v[i+2][j]+8.0*(v[i+1][j]-v[i-1][j])+v[i-2][j])/(12.0*dx) 
                         -abs(u[i][j])*(v[i+2][j]-4.0*v[i+1][j]+6.0*v[i][j]-4.0*v[i-1][j]+v[i-2][j])/(4.0*dx);
        }
    }
    // advection term in y-direction(y方向についても同様の計算を行う)
    for (int i=i_1+1; i<=i_2-1; i++){
        u[i][j_2+1] = 2.0*u[i][j_2]-u[i][j_2-1];
        u[i][j_3-1] = 2.0*u[i][j_3]-u[i][j_3+1];
        v[i][j_2+1] = 2.0*v[i][j_2]-v[i][j_2-1];
        v[i][j_3-1] = 2.0*v[i][j_3]-v[i][j_3+1];
    }

    for (int i=i_2+1; i<=i_3-1; i++){
        u[i][j_1+1] = 2.0*u[i][j_1]-u[i][j_1-1];
        u[i][j_4-1] = 2.0*u[i][j_4]-u[i][j_4+1];
        v[i][j_1+1] = 2.0*v[i][j_1]-v[i][j_1-1];
        v[i][j_4-1] = 2.0*v[i][j_4]-v[i][j_4+1];
    }

    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            urhs[i][j] = urhs[i][j]-v[i][j]*(-u[i][j+2]+8.0*(u[i][j+1]-u[i][j-1])+u[i][j-2])/(12.0*dy)
                         -abs(v[i][j])*(u[i][j+2]-4.0*u[i][j+1]+6.0*u[i][j]-4.0*u[i][j-1]+u[i][j-2])/(4.0*dy);
			vrhs[i][j] = vrhs[i][j]-v[i][j]*(-v[i][j+2]+8.0*(v[i][j+1]-v[i][j-1])+v[i][j-2])/(12.0*dy) 
                         -abs(v[i][j])*(v[i][j+2]-4.0*v[i][j+1]+6.0*v[i][j]-4.0*v[i][j-1]+v[i][j-2])/(4.0*dy);
        }
    }
    // update
    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            u[i][j] = u[i][j]+dt*urhs[i][j];
            v[i][j] = v[i][j]+dt*vrhs[i][j];
        }
    }
}

void calc_omega(vd2 & u, vd2 & v){
    for (int i=2; i<=mx-1; i++){
        for (int j=2; j<=my-1; j++){
            // 物体内部は計算しない
            if (i>i_1 && i<=i_2 && j>j_2 && j<j_3){
                continue;
            }
            else if (i>i_2 && i<i_3 && j>j_1 && j<j_4){
                continue;
            }
            double uy = (u[i][j+1]-u[i][j-1])/(2.0*dy);
            double vx = (v[i+1][j]-v[i-1][j])/(2.0*dx);
            omega[i][j] = vx-uy;
        }
    }
}

// 流れをステップごとに時間発展させ，ファイルに保存する
void solve_flow(vd1 & x, vd1 & y, vd2 & u, vd2 & v, vd2 & p){
    // 書き込みファイルの作成
    std::ofstream MAC_deco_1;
    MAC_deco_1.open("MAC_deco_1000.csv", std::ios::trunc);
    std::ofstream MAC_deco_2;
    MAC_deco_2.open("MAC_deco_2000.csv", std::ios::trunc);
    std::ofstream MAC_deco_3;
    MAC_deco_3.open("MAC_deco_3000.csv", std::ios::trunc);
    std::ofstream MAC_deco_4;
    MAC_deco_4.open("MAC_deco_4000.csv", std::ios::trunc);
    std::ofstream MAC_deco_5;
    MAC_deco_5.open("MAC_deco_5000.csv", std::ios::trunc);
    std::ofstream MAC_deco_omega;
    MAC_deco_omega.open("MAC_deco_omega.csv", std::ios::trunc);
    std::ofstream rhs_hs;
    rhs_hs.open("history_of_rhs.csv", std::ios::trunc);
    std::ofstream rhs_hs2;
    rhs_hs2.open("history_of_rhs2.csv", std::ios::trunc);
    std::ofstream dimless_coef_hs;
    dimless_coef_hs.open("history_dimless_coef.csv", std::ios::trunc);
    //  時間計測用変数を確保
    std::chrono::system_clock::time_point start_time, end_time; 
    // 計測スタート
    start_time = std::chrono::system_clock::now();
    for (int n=1; n<=nlast; n++){
        simu_time = simu_time+dt;
        poiseq(u, v, p);
        bcforp(p);
        veloeq(u, v, p);
        bcforv(u, v);
        calc_omega(u, v);

        double cd = 0.0;
        for (int j=j_2; j<=j_3-1; j++){
            double cpfore1 = (2.0*p[i_1][j]+2.0*p[i_1][j+1])/2.0;
            // double cpback = (2.0*p[i_2][j]+2.0*p[i_2][j+1])/2.0;
            cd = cd+cpfore1*dy;
        }
        for (int j=j_1; j<=j_2-1; j++){
            double cpfore2 = (2.0*p[i_2][j]+2.0*p[i_2][j+1])/2.0;
            cd = cd+cpfore2*dy;
        }
        for (int j=j_3; j<=j_4-1; j++){
            double cpfore3 = (2.0*p[i_2][j]+2.0*p[i_2][j+1])/2.0;
            cd = cd+cpfore3*dy;
        }
        for (int j=j_1; j<=j_4-1; j++){
            double cpback = (2.0*p[i_3][j]+2.0*p[i_3][j+1])/2.0;
            cd = cd-cpback*dy;
        }

        double cl = 0.0;
        for (int i=i_1; i<=i_2-1; i++){
            double cpbtm1 = (2.0*p[i][j_2]+2.0*p[i+1][j_2])/2.0;
            double cptop1 = (2.0*p[i][j_3]+2.0*p[i+1][j_3])/2.0;
            cl = cl+(cpbtm1-cptop1)*dx;
        }
        for (int i=i_2; i<=i_3-1; i++){
            double cpbtm2 = (2.0*p[i][j_1]+2.0*p[i+1][j_1])/2.0;
            double cptop2 = (2.0*p[i][j_4]+2.0*p[i+1][j_4])/2.0;
            cl = cl+(cpbtm2-cptop2)*dx;
        }

        double cp1 = 2.0*p[i_3+i_2-i_1][j_1];
        double cp2 = 2.0*p[i_3+i_2-i_1][j_4];

        end_time = std::chrono::system_clock::now();
        double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count()/1000.0; 
        std::cout << "step : " << n <<  "/" << nlast << " elapsed_time : " << elapsed_time << "[s]" << std::endl;
        std::cout << "elapsed time in simulation :" << dt*n << " [s]" << std::endl;
        // ステップごとのrhsと各抵抗係数を書き込み
        rhs_hs << simu_time << "," << rhs[int((i_1+i_2)/2.0)+50][int((j_1+j_2)/2.0)] << std::endl;
        rhs_hs2 << simu_time << "," << rhs[int((i_1+i_2)/2.0)+150][int((j_1+j_2)/2.0)] << std::endl;
        dimless_coef_hs << simu_time << "," << cd << "," << cl << "," << cp1 << "," << cp2 << std::endl;
        // csvに，代表的なステップにおける圧力分布を書き込み
        if (n == 1000){
            for (int i=1; i<=mx; i++){
                for (int j=1; j<=my; j++){
                    MAC_deco_1 << x[i] << "," << y[j] << "," << 2.0*p[i][j] << std::endl; // pの値を2倍することに注意
                }
                MAC_deco_1 << std::endl; // gnuplotのpm3dのために必要
            }
        }

        if (n == 2000){
            for (int i=1; i<=mx; i++){
                for (int j=1; j<=my; j++){
                    MAC_deco_2 << x[i] << "," << y[j] << "," << 2.0*p[i][j] << std::endl; // pの値を2倍することに注意
                }
                MAC_deco_2 << std::endl; // gnuplotのpm3dのために必要
            }
        }

        if (n == 3000){
            for (int i=1; i<=mx; i++){
                for (int j=1; j<=my; j++){
                    MAC_deco_3 << x[i] << "," << y[j] << "," << 2.0*p[i][j] << std::endl; // pの値を2倍することに注意
                }
                MAC_deco_3 << std::endl; // gnuplotのpm3dのために必要
            }
        }

        if (n == 4000){
            for (int i=1; i<=mx; i++){
                for (int j=1; j<=my; j++){
                    MAC_deco_4 << x[i] << "," << y[j] << "," << 2.0*p[i][j] << std::endl; // pの値を2倍することに注意
                }
                MAC_deco_4 << std::endl; // gnuplotのpm3dのために必要
            }
        }

        if (n == nlast){
            for (int i=1; i<=mx; i++){
                for (int j=1; j<=my; j++){
                    MAC_deco_5 << x[i] << "," << y[j] << "," << 2.0*p[i][j] << std::endl; // pの値を2倍することに注意
                    MAC_deco_omega << x[i] << "," << y[j] << "," << omega[i][j] << std::endl;
                }
                MAC_deco_5 << std::endl; // gnuplotのpm3dのために必要
                MAC_deco_omega << std::endl; // gnuplotのpm3dのために必要
            }
        }
    }
    std::cout << "simu_time[s] : " << simu_time << std::endl;
    MAC_deco_1.close();
    MAC_deco_2.close();
    MAC_deco_3.close();
    MAC_deco_4.close();
    MAC_deco_5.close();
    MAC_deco_omega.close();
    dimless_coef_hs.close();
    rhs_hs.close();
    rhs_hs2.close();
}