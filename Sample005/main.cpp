/*
 * convection eq.
 *
 * CIP method.
 * 
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>

constexpr int N = 64;
constexpr double Lx =  2.0;
constexpr double Ly =  2.0;
constexpr double xstart = -1.0;
constexpr double ystart = -1.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.0001;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 2.0 * pi;
constexpr int INTV = 100;
constexpr int plus = 2;

namespace mino2357{
    template <typename T = double>
    class extendedArray{
    private:
        T* u;
        int num;
    public:
        extendedArray(int n){
            num = n;
            u = new T[num * num];
        }
        ~extendedArray(){ delete[] u;}

        constexpr T& operator()(int i, int j){
            return u[(((j % num) + num) % num) * num + ((i % num) + num) % num];
        }
    };
    
    template <typename T = double>
    T x(int i){
        return i;
        return xstart + i * dx;
    }
    
    template <typename T = double>
    T y(int j){
        return j;
        return ystart + j * dy;
    }

    template <typename T = double>
    T f(T x, T y){
        return -y;
        return std::sin(pi * (x + y)) + 2.0;
    }
    
    template <typename T = double>
    T g(T x, T y){
        return x;
        return std::sin(pi * (x + y)) + 2.0;
    }

    template <typename T = double>
    void makeInitFuncU(extendedArray<T>& u){
       // T a = 0.5;
       // T b = 0.5;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(xstart+i*dx >= 0.2 && xstart+i*dx <= 0.6 && ystart+j*dy >= 0.2 && ystart+j*dy <= 0.6){
                    u(i, j) = 1.0;
                }else{
                    u(i, j) = 0.0;
                }
                //u(i, j) = //std::sin(pi * (xstart + i * dx + ystart + j * dy)) + 2.0;
                  //        std::exp(- 25.0 * (std::pow(xstart + i * dx - a, 2) + std::pow(ystart + j * dy - b, 2)));
            }
        }
    }

    template <typename T = double>
    void makeInitFuncV(extendedArray<T>& u){
        T a = 0.5;
        T b = 0.5;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                /*
                if(xstart+i*dx >= 0.2 && xstart+i*dx <= 0.6 && ystart+j*dy >= 0.2 && ystart+j*dy <= 0.6){
                    u(i, j) = 1.0;
                }else{
                    u(i, j) = 0.0;
                }
                */
                u(i, j) = //std::sin(pi * (xstart + i * dx + ystart + j * dy)) + 2.0;
                          std::exp(- 25.0 * (std::pow(xstart + i * dx - a, 2) + std::pow(ystart + j * dy - b, 2)));
            }
        }
    }

    template <typename T = double>
    void makeInitGradX(extendedArray<T>& u, extendedArray<T>& gx1){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gx1(i, j) = (u(i+1, j) - u(i-1, j)) / (2.0 * dx);
            }
        }
    }

    template <typename T = double>
    void makeInitGradY(extendedArray<T>& u, extendedArray<T>& gy1){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gy1(i, j) = (u(i, j+1) - u(i, j-1)) / (2.0 * dy);
            }
        }
    }

    namespace coeff{
        namespace X{
            namespace U{
                template< typename T = double>
                T a(int i, int j, extendedArray<T>& u, extendedArray<T>& ugx){
                    if(u(x(i), y(j)) >= 0.0){
                        return (ugx(i, j) + ugx(i-1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i-1, j)) / (dx * dx * dx);
                    }
                    return (ugx(i, j) + ugx(i+1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i+1, j)) / (dx * dx * dx);
                }
            
                template< typename T = double>
                T b(int i, int j, extendedArray<T>& u, extendedArray<T>& ugx){
                    if(u(x(i), y(j)) >= 0.0){
                        return 3.0 * (u(i-1, j) - u(i, j)) / (dx * dx) + (2.0 * ugx(i, j) + ugx(i-1, j)) / (dx);
                    }
                    return 3.0 * (u(i+1, j) - u(i, j)) / (dx * dx) + (2.0 * ugx(i, j) + ugx(i+1, j)) / (dx);
                }
            
                template< typename T = double>
                T c(int i, int j, extendedArray<T>& g){
                    return g(i, j);
                }
            
                template< typename T = double>
                T d(int i, int j, extendedArray<T>& u){
                    return u(i, j);
                }
            }

            namespace V{
                template< typename T = double>
                T a(int i, int j, extendedArray<T>& u, extendedArray<T>& v, extendedArray<T>& ugx){
                    if(v(x(i), y(j)) >= 0.0){
                        return (ugx(i, j) + ugx(i-1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i-1, j)) / (dx * dx * dx);
                    }
                    return (ugx(i, j) + ugx(i+1, j)) / (dx * dx) - 2.0 * (u(i, j) - u(i+1, j)) / (dx * dx * dx);
                }
            
                template< typename T = double>
                T b(int i, int j, extendedArray<T>& u, extendedArray<T>& v, extendedArray<T>& ugx){
                    if(v(x(i), y(j)) >= 0.0){
                        return 3.0 * (u(i-1, j) - u(i, j)) / (dx * dx) + (2.0 * ugx(i, j) + ugx(i-1, j)) / (dx);
                    }
                    return 3.0 * (u(i+1, j) - u(i, j)) / (dx * dx) + (2.0 * ugx(i, j) + ugx(i+1, j)) / (dx);
                }
            
                template< typename T = double>
                T c(int i, int j, extendedArray<T>& g){
                    return g(i, j);
                }
            
                template< typename T = double>
                T d(int i, int j, extendedArray<T>& u){
                    return u(i, j);
                }
            }
        }

        namespace Y{
            namespace U{
                template< typename T = double>
                T a(int i, int j, extendedArray<T>& u, extendedArray<T>& gy){
                    if(u(x(i), y(j)) >= 0.0){
                        return (gy(i, j) + gy(i, j-1)) / (dy * dy) - 2.0 * (u(i, j) - u(i, j-1)) / (dy * dy * dy);
                    }
                    return (gy(i, j) + gy(i, j+1)) / (dy * dy) - 2.0 * (u(i, j) - u(i, j+1)) / (dy * dy * dy);
                }
            
                template< typename T = double>
                T b(int i, int j, extendedArray<T>& u, extendedArray<T>& gy){
                    if(u(x(i), y(j)) >= 0.0){
                        return 3.0 * (u(i, j-1) - u(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j-1)) / (dy);
                    }
                    return 3.0 * (u(i, j+1) - u(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j+1)) / (dy);
                }
            
                template< typename T = double>
                T c(int i, int j, extendedArray<T>& g){
                    return g(i, j);
                }
            
                template< typename T = double>
                T d(int i, int j, extendedArray<T>& u){
                    return u(i, j);
                }
            }

            namespace V{
                template< typename T = double>
                T a(int i, int j, extendedArray<T>& v, extendedArray<T>& gy){
                    if(v(x(i), y(j)) >= 0.0){
                        return (gy(i, j) + gy(i, j-1)) / (dy * dy) - 2.0 * (v(i, j) - v(i, j-1)) / (dy * dy * dy);
                    }
                    return (gy(i, j) + gy(i, j+1)) / (dy * dy) - 2.0 * (v(i, j) - v(i, j+1)) / (dy * dy * dy);
                }
            
                template< typename T = double>
                T b(int i, int j, extendedArray<T>& v, extendedArray<T>& gy){
                    if(v(x(i), y(j)) >= 0.0){
                        return 3.0 * (v(i, j-1) - v(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j-1)) / (dy);
                    }
                    return 3.0 * (v(i, j+1) - v(i, j)) / (dy * dy) + (2.0 * gy(i, j) + gy(i, j+1)) / (dy);
                }
            
                template< typename T = double>
                T c(int i, int j, extendedArray<T>& g){
                    return g(i, j);
                }
            
                template< typename T = double>
                T d(int i, int j, extendedArray<T>& u){
                    return u(i, j);
                }
            }
        }
    }

    template <typename T = double>
    void succUUU(extendedArray<T>& u1, extendedArray<T>& ugx, extendedArray<T>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(u1(x(i), y(j)) >= 0){
                    up = - u1(x(i), y(j)) * dt;
                }else{
                    up = u1(x(i), y(j)) * dt;
                }
                u2(i, j) =   coeff::X::U::a(i, j, u1, ugx) * up * up * up 
                           + coeff::X::U::b(i, j, u1, ugx) * up * up
                           + coeff::X::U::c(i, j, ugx)     * up
                           + coeff::X::U::d(i, j, u1);
            }
        }
    }

    template <typename T = double>
    void succUVU(extendedArray<T>& u1, extendedArray<T>& v1, extendedArray<T>& ugy, extendedArray<T>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(v1(x(i), y(j)) >= 0){
                    up = - v1(x(i), y(j)) * dt;
                }else{
                    up = v1(x(i), y(j)) * dt;
                }
                u2(i, j) =   coeff::Y::V::a(i, j, u1, v1, ugy) * up * up * up 
                           + coeff::Y::V::b(i, j, u1, v1, ugy) * up * up
                           + coeff::Y::V::c(i, j, ugy)     * up
                           + coeff::Y::V::d(i, j, u1);
            }
        }
    }

    template <typename T = double>
    void succVUV(extendedArray<T>& v1, extendedArray<T>& u1, extendedArray<T>& vgx, extendedArray<T>& v2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(u1(x(i), y(j)) >= 0){
                    up = - v1(x(i), y(j)) * dt;
                }else{
                    up = v1(x(i), y(j)) * dt;
                }
                v2(i, j) =   coeff::X::U::a(i, j, u1, vgx) * up * up * up 
                           + coeff::X::U::b(i, j, u1, vgx) * up * up
                           + coeff::X::U::c(i, j, vgx)     * up
                           + coeff::X::U::d(i, j, u1);
            }
        }
    }
    
    template <typename T = double>
    void succVVV(extendedArray<T>& v1, extendedArray<T>& vgy, extendedArray<T>& v2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(v1(x(i), y(j)) >= 0){
                    up = - v1(x(i), y(j)) * dt;
                }else{
                    up = v1(x(i), y(j)) * dt;
                }
                v2(i, j) =   coeff::Y::V::a(i, j, v1, vgy) * up * up * up 
                           + coeff::Y::V::b(i, j, v1, vgy) * up * up
                           + coeff::Y::V::c(i, j, vgy)     * up
                           + coeff::Y::V::d(i, j, v1);
            }
        }
    }
/*
    template <typename T = double>
    void succGX(extendedArray<T>& u, extendedArray<T>& gx1, extendedArray<T>& gx2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(f(x(i), y(j)) >= 0){
                    up = - f(x(i), y(j)) * dt;
                }else{
                    up = f(x(i), y(j)) * dt;
                }
                gx2(i, j) =   3.0 * coeff::X::a(i, j, u, gx1) * up * up
                            + 2.0 * coeff::X::b(i, j, u, gx1) * up
                            +       coeff::X::c(i, j, gx1);
            }
        }
    }
    
    template <typename T = double>
    void succGY(extendedArray<T>& u, extendedArray<T>& gy1, extendedArray<T>& gy2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(g(x(i), y(j)) >= 0){
                    up = - g(x(i), y(j)) * dt;
                }else{
                    up = g(x(i), y(j)) * dt;
                }
                gy2(i, j) =   3.0 * coeff::Y::a(i, j, u, gy1) * up * up
                            + 2.0 * coeff::Y::b(i, j, u, gy1) * up
                            +       coeff::Y::c(i, j, gy1);
            }
        }
    }   */
}


int main(){
    auto u1   = mino2357::extendedArray<>(N + 1);
    auto u2   = mino2357::extendedArray<>(N + 1);
    auto u3   = mino2357::extendedArray<>(N + 1);
    auto ugx1 = mino2357::extendedArray<>(N + 1);
    auto ugy1 = mino2357::extendedArray<>(N + 1);
    auto ugx2 = mino2357::extendedArray<>(N + 1);
    auto ugy2 = mino2357::extendedArray<>(N + 1);
    auto ugx3 = mino2357::extendedArray<>(N + 1);
    auto ugy3 = mino2357::extendedArray<>(N + 1);

    auto v1   = mino2357::extendedArray<>(N + 1);
    auto v2   = mino2357::extendedArray<>(N + 1);
    auto v3   = mino2357::extendedArray<>(N + 1);
    auto vgx1 = mino2357::extendedArray<>(N + 1);
    auto vgy1 = mino2357::extendedArray<>(N + 1);
    auto vgx2 = mino2357::extendedArray<>(N + 1);
    auto vgy2 = mino2357::extendedArray<>(N + 1);
    auto vgx3 = mino2357::extendedArray<>(N + 1);
    auto vgy3 = mino2357::extendedArray<>(N + 1);

    mino2357::makeInitFuncU(u1);
    mino2357::makeInitFuncV(v1);
    
    mino2357::makeInitGradX(u1, ugx1);
    mino2357::makeInitGradY(u1, ugy1);
    
    mino2357::makeInitGradX(v1, vgx1);
    mino2357::makeInitGradY(v1, vgy1);

    std::cout << dt / dx << std::endl;

    /***********************************************/
    FILE* gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set pm3d\n");
    //fprintf(gp, "set pm3d map\n");
    fprintf(gp, "set contour\n");
    fprintf(gp, "set xr [%f:%f]\n", xstart, xstart + Lx);
    fprintf(gp, "set yr [%f:%f]\n", ystart, ystart + Ly);
  	//fprintf(gp, "set zr [-0.1:1.1]\n");
    fprintf(gp, "set size square\n");

    fprintf(gp, "splot '-' w l\n");
    for(int i=0; i<=N; i = i + plus){
        for(int j=0; j<=N; j = j + plus){
            fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1(i, j));
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
    fflush(gp);
    /******************************************************/

    std::cout << "Enter!" << std::endl;
    getchar();

    double t = 0.0;

    for(int it=0; t<=tLimit; ++it){
        
        mino2357::succUUU(u1, ugx1, u2);
        mino2357::succUVU(u2, v1, ugy1, u3);
        mino2357::succVUV(v1, u1, vgx1, v2);
        mino2357::succVVV(v2, vgy1, v3);
        /*
        mino2357::succGX(u1, ugx1, ugx2);
        mino2357::succGY(u1, ugx1, ugx3);
        mino2357::succGX(u1, ugy1, ugy2);
        mino2357::succGY(u1, ugy1, ugy3);
        */
		for(int i=0; i<=N; ++i){
			for(int j=0; j<=N; ++j){
				ugx3(i, j) = 0.5 * (ugx2(i, j) + ugx3(i, j));
				ugy3(i, j) = 0.5 * (ugy2(i, j) + ugy3(i, j));
			}
		}

        /*********************************************/
        if(it%INTV == 0){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1(i, j));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        /*********************************************/

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u1(i, j)  = u3(i, j);
                ugx1(i, j) = ugx3(i, j);
                ugy1(i, j) = ugy3(i, j);
            }
        }

        t += dt;
    }
    fprintf(gp, "splot '-' w l\n");
    for(int i=0; i<=N; i = i + plus){
      	for(int j=0; j<=N; j = j + plus){
           	fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1(i, j));
      	}
       	fprintf(gp, "\n");
 	}
   	fprintf(gp, "e\n");
   	fflush(gp);

	std::cout << "simulation done.(mean)" << std::endl;
	getchar();

    pclose(gp);
}
