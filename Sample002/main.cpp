/*
 * 拡散方程式
 *
 * Takaaki MINOMO
 *
 */

#include <iostream>
#include <array>
#include <cmath>
#include <sprout/cmath.hpp> 

constexpr int N = 128;
constexpr double Lx = 1.0;
constexpr double Ly = 1.0;
constexpr double cx = 1.0;
constexpr double cy = 1.0;
constexpr double D = 0.0;
constexpr double tLimit = 1000;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.002;
constexpr double x  = 0.;
constexpr double y  = 0.;
constexpr int INTV = 3;

namespace mino2357{
    
    template <typename T = double>
    class extendedArray{
    private:
        T* u;
        int Num;
    public:
        extendedArray(int n){ Num = n; u = new T[(Num + 1) * (Num + 1)];}
        ~extendedArray(){ delete[] u;}
        
        constexpr T& operator()(int, int);

    };
    
    template <typename T>
    constexpr T& extendedArray<T>::operator()(int i, int j){
        return u[(j%Num) * Num + (i%Num)];
    }

    template <typename T = double>
    constexpr T initFunc(T x, T y) noexcept {
        T a = 5.0;
        T b = 5.0;
        return std::exp(- 25.0 *((x - a * Lx / 10.0) * (x - a * Lx / 10.0) + (y - b * Ly / 10.0) * (y - b * Ly / 10.0)));
    }

    template <typename T = double>
    constexpr void makeInitFunc(extendedArray<T>& u) noexcept {
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u(i, j) = initFunc<T>(i * dx, j * dx);
            }
        }
    }
}

int main(){
    double t = 0.0;	

    auto u1 = mino2357::extendedArray<>(N);
    auto u2 = mino2357::extendedArray<>(N);
    mino2357::makeInitFunc<>(u1);
  
    std::cout << D * dt / (dx * dx) << std::endl;
    if(cx * dt / dx >=1.0 || cy * dt / dy >= 1.0 || D * dt / (dx * dx) > 0.5 || D * dt / (dy * dy) > 0.5){
        std::cout << "安定性条件を満たしていません．" << std::endl;
    }

    /**********************************************************************/
    /*                 可視化の設定(gnuplot)                              */
    /**********************************************************************/
    
    std::FILE *gp = popen( "gnuplot -persist", "w" );
    fprintf(gp, "set pm3d\n");
    //fprintf(gp, "set pm3d map\n");
    fprintf(gp, "set contour\n");
    fprintf(gp, "set xr [0:%f]\n", Lx);
    fprintf(gp, "set yr [0:%f]\n", Ly);
    fprintf(gp, "set zr [0.0:1.0]\n");
    fprintf(gp, "set size square\n");
    //fprintf(gp, "set grid\n");
    //fprintf(gp, "unset key\n");
    
    //初期条件描画
    fprintf(gp, "splot '-'w l\n");
    for(int i=0; i<=N; ++i){
        for(int j=0; j<=N; ++j){
            fprintf(gp, "%f %f %f\n", x + i * dx, y + j * dy, u1(i, j));
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fflush(gp);

    std::cout << "Enterキーを押してください．" << std::endl;
    getchar();

    //タイムループ
    for(int it = 0; t<tLimit; ++it) {   

        //拡散
        /*
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u2(i, j) = u1(i, j) + D * dt / (dx * dx) * (u1(i-1, j) - 2.0 * u1(i, j) + u1(i+1, j)) + D * dt / (dy * dy) * (u1(i, j-1) - 2.0 * u1(i, j) + u1(i, j+1)); 
            }        
        }
        */

        //移流 cx > 0 cy > 0 の対応しかまだしてない．
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u2(i, j) = u1(i,j) - cx * dt / dx * (u1(i,j) - u1(i-1,j)) - cy * dt / dy * (u1(i, j) - u1(i, j-1));
            }
        }
        
        //u2の描画
        if(it%INTV == 0){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; ++i){
                for(int j=0; j<=N; ++j){
                    fprintf(gp, "%f %f %f\n", x + i * dx, y + j * dy, u2(i, j));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
    
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u1(i, j) = u2(i, j);
            }        
        }

    }

    //FILEポインタの解放
    pclose(gp);
}
