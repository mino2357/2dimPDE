/*
 * 非線形移流方程式
 *
 * M-CIP method.
 *
 * Takaaki MINOMO
 *
 */

#include <iostream>
#include <cmath>

//parameter
constexpr int N = 128;
constexpr double Lx = 1.0;
constexpr double Ly = 1.0;
constexpr double cx = 1.0;
constexpr double cy = 1.0;
constexpr double D = 0.0;
constexpr double tLimit = 1000;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.001;
constexpr double x  = 0.0;
constexpr double y  = 0.0;
constexpr double pi = 3.14159265358979323;
constexpr int INTV = 10;
constexpr int plus = 4;

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
        return std::sin(2.0 * pi * x);// + std::cos(2.0 * pi * y);
        //if(x >= 0.3 && x <= 0.8 && y >= 0.3 && y <= 0.8){
        if(x >= 0.5){
            return 0.0;
        }
        return 1.0;
        
        T a = 5.0;
        T b = 5.0;
        return std::exp( - 25.0 *((x - a / 10.0) * (x - a / 10.0)
                                + (y - b / 10.0) * (y - b / 10.0)));
    }

    template <typename T = double>
    constexpr void makeInitFunc(extendedArray<T>& u) noexcept {
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u(i, j) = initFunc<T>(x + i * dx, y + j * dx);
            }
        }
    }

    template <typename T = double>
    constexpr void makeInitGradX(extendedArray<T>& u, extendedArray<T>& gx){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                //g(i, j) = (u(i-1, j) - u(i+1, j))/ (2.0 * dx) + (u(i, j-1) - u(i, j+1)) / (2.0 * dy);
                gx(i, j) = (u(i-1, j) - u(i+1, j))/ (2.0 * dx);
            }
        }
    }
    
    template <typename T = double>
    constexpr void makeInitGradY(extendedArray<T>& u, extendedArray<T>& gy){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gy(i, j) = (u(i, j-1) - u(i, j+1)) / (2.0 * dy);
            }
        }
    }
    
    /******************************************************************************************/

    namespace coeff{
        namespace X{
            template <typename T = double>
            T a(int i, int j, extendedArray<T>& u, extendedArray<T>& gx){
                if(u(i, j) * cx >= 0.0){
                    return (gx(i, j) + gx(i-1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i-1, j)) / (dx * dx * dx);
                }else{
                    return (gx(i, j) + gx(i+1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i+1, j)) / (dx * dx * dx);
                }
            }

            template <typename T = double>
            T b(int i, int j, extendedArray<T>& u, extendedArray<T>& gx){
                if(cx * u(i, j) >= 0.0){
                    return 3.0 * (u(i-1, j) - u(i, j)) / (dx * dx) + (2.0 * gx(i, j) + gx(i-1, j)) / (dx);
                }else{
                    return 3.0 * (u(i+1, j) - u(i, j)) / (dx * dx) + (2.0 * gx(i, j) + gx(i+1, j)) / (dx);
                }
            }

            template <typename T = double>
            T c(int i, int j, extendedArray<T>& gx){
                return gx(i, j);
            }

            template <typename T = double>
            T d(int i, int j, extendedArray<T>& u){
                return u(i, j);
            }
        }
        namespace Y{
            template <typename T = double>
            T a(int i, int j, extendedArray<T>& u, extendedArray<T>& gy){
                if(cy * u(i, j) >= 0.0){
                    return (gy(i, j) + gy(i, j-1)) / (dy * dy) - 2.0 * (u(i, j) - u(i, j-1)) / (dy * dy * dy);
                }else{
                    return (gy(i, j) + gy(i, j+1)) / (dy * dy) - 2.0 * (u(i, j) - u(i, j+1)) / (dy * dy * dy);
                }
            }

            template <typename T = double>
            T b(int i, int j, extendedArray<T>& u, extendedArray<T>& gy){
                if(cy * u(i, j) >= 0.0){
                    return 3.0 * (u(i, j-1) - u(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j-1)) / (dy);
                }else{
                    return 3.0 * (u(i, j+1) - u(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j+1)) / (dy);
                }
            }

            template <typename T = double>
            T c(int i, int j, extendedArray<T>& gy){
                return gy(i, j);
            }

            template <typename T = double>
            T d(int i, int j, extendedArray<T>& u){
                return u(i, j);
            }
        }
    }
    
    template <typename T = double>
    void makeSuccUX(extendedArray<T>& u, extendedArray<T>& gx, extendedArray<T>& succUX){
        T zx;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(cx * u(i, j) >= 0.0){
                    zx = - u(i, j) * cx * dt;
                }else{
                    zx = u(i, j) * cx * dt;
                }
                succUX(i, j) = coeff::X::a<>(i, j, u, gx) * zx * zx * zx + coeff::X::b<>(i, j, u, gx) * zx * zx + coeff::X::c<>(i, j, gx) * zx + coeff::X::d<>(i, j, u);
            }
        }
    }

    template <typename T = double>
    void makeSuccUY(extendedArray<T>& u, extendedArray<T>& gy, extendedArray<T>& succUY){
        T zy;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(cy * u(i, j) >= 0.0){
                    zy = - u(i, j) * cy * dt;
                }else{
                    zy = u(i, j) * cy * dt;
                }
                succUY(i, j) = coeff::Y::a<>(i, j, u, gy) * zy * zy * zy + coeff::Y::b<>(i, j, u, gy) * zy * zy + coeff::Y::c<>(i, j, gy) * zy + coeff::Y::d<>(i, j, u);
            }
        }
    }

    template <typename T = double>
    void makeSuccGX(extendedArray<T>& u, extendedArray<T>& gx, extendedArray<T>& succGX){
        T zx;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(cx * u(i, j) >= 0.0){
                    zx = - u(i,j) * cx * dt;
                }else{
                    zx = u(i,j) * cx * dt;
                }
                succGX(i, j) = 3.0 * coeff::X::a(i, j, u, gx) * zx * zx + 2.0 * coeff::X::b(i, j, u, gx) * zx + coeff::X::c(i, j, gx);
            }
        }
    }

    template <typename T = double>
    void makeSuccGY(extendedArray<T>& u, extendedArray<T>& gy, extendedArray<T>& succGY){
        T zy;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(cy * u(i, j) >= 0.0){
                    zy = - u(i,j) * cy * dt;
                }else{
                    zy = u(i,j) * cy * dt;
                }
                succGY(i, j) = 3.0 * coeff::Y::a(i, j, u, gy) * zy * zy + 2.0 * coeff::Y::b(i, j, u, gy) * zy + coeff::Y::c(i, j, gy);
            }
        }
    }
}

int main(){
    double t = 0.0;	

    auto u1 = mino2357::extendedArray<>(N);
    auto ux = mino2357::extendedArray<>(N);
    auto uy = mino2357::extendedArray<>(N);
    auto u2 = mino2357::extendedArray<>(N);
    auto gx1 = mino2357::extendedArray<>(N);
    auto gx2 = mino2357::extendedArray<>(N);
    auto gy1 = mino2357::extendedArray<>(N);
    auto gy2 = mino2357::extendedArray<>(N);
    auto gx = mino2357::extendedArray<>(N);
    auto gy = mino2357::extendedArray<>(N);
    
    mino2357::makeInitFunc<>(u1);
    mino2357::makeInitGradX<>(u1, gx1);
    mino2357::makeInitGradY<>(u1, gy1);
  
    std::cout << cx * dt / dx << std::endl;
    if(cx * dt / dx >= 1.0 || cy * dt / dy >= 1.0 || D * dt / (dx * dx) > 0.5 || D * dt / (dy * dy) > 0.5){
        std::cout << "Not satisfy the stability condition!" << std::endl;
    }

    /**********************************************************************/
    /*                visualization(gnuplot)                              */
    /**********************************************************************/
    
    std::FILE *gp = popen( "gnuplot -persist", "w" );
    fprintf(gp, "set pm3d\n");
    //fprintf(gp, "set pm3d map\n");
    fprintf(gp, "set contour\n");
    fprintf(gp, "set xr [%f:%f]\n", x, x + Lx);
    fprintf(gp, "set yr [%f:%f]\n", y, y + Ly);
    //fprintf(gp, "set zr [-0.1:1.1]\n");
    fprintf(gp, "set size square\n");
    //fprintf(gp, "set grid\n");
    //fprintf(gp, "unset key\n");
    
    //first condition visualization
    fprintf(gp, "splot '-'w l\n");
    for(int i=0; i<=N; i=i+plus){
        for(int j=0; j<=N; j=j+plus){
            fprintf(gp, "%f %f %f\n", x + i * dx, y + j * dy, u1(i, j));
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fflush(gp);

    std::cout << "Enter!" << std::endl;
    getchar();

    //Time Loop
    for(int it = 0; t<tLimit; ++it) {   

        //Calc next u and g
        mino2357::makeSuccUX(u1, gx1, ux);
        mino2357::makeSuccUY(u1, gy1, uy);

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u2(i, j) = ux(i, j) + 0.5 * uy(i, j);
            }
        }

        mino2357::makeSuccGX(u1, gx1, gx);
        mino2357::makeSuccGY(u1, gy1, gy);
       
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gx2(i, j) = gx(i, j) - cx * dt * gx(i, j) * gx(i, j);
                gy2(i, j) = gy(i, j) - cy * dt * gy(i, j) * gy(i, j);
            }
        }

        //u2: visualization
        if(it%INTV == 0){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i=i+plus){
                for(int j=0; j<=N; j=j+plus){
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
                gx1(i, j) = gx2(i, j);
                gy1(i, j) = gy2(i, j);
            }        
        }

    }

    //free file pointer
    pclose(gp);
}
