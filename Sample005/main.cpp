/*
 * NS eq.
 *
 * CIP method.
 * 
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>
#include <array>

constexpr int N = 256;
constexpr double Re = 100.0;
constexpr double Lx =  2.0;
constexpr double Ly =  2.0;
constexpr double xstart = - 1.0;
constexpr double ystart = - 1.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.00001;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 20000.0 * pi;
constexpr int INTV = 100;
constexpr int plus = 8;

namespace rittai3d{
	namespace utility{
		// [ minimum, maximum ) の範囲でラップアラウンド
		template <typename T>
		constexpr T wrap_around( T value, T minimum, T maximum )
		{
			const T n = ( value - minimum ) % ( maximum - minimum );
			return n >= 0 ? ( n + minimum ) : ( n + maximum ); 
		}
	}
}

namespace sksat {
	template<std::size_t NN, typename T = double>
	class array_wrapper {
	private:
		std::array<T,NN> arr;
	public:
		constexpr array_wrapper() : arr() {}
		~array_wrapper() = default;

		T& operator[](int i){
			i = rittai3d::utility::wrap_around(i,0,static_cast<int>(NN));
			return arr[i];
		}
	};
}

namespace mino2357{
    
    template <typename T = double>
    void makeInitFuncU(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u){
        T x, y;
        T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.0 && x <= 0.4 && y >= 0.0 && y <= 0.4){
                    u[i][j] = 1.0;
                }else{
                    u[i][j] = 0.0;
                }
                
               u[i][j] = 0.1 * (std::sin(pi * (x)));// + std::sin(pi * (y));
                       //0.1 * std::exp(- 20.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
                
            } 
        }
    }

    template <typename T = double>
    void makeInitFuncV(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& v){
        
        T x, y;
        T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
              
                if(x >= 0.0 && x <= 0.4 && y >= 0.0 && y <= 0.4){
                    v[i][j] = 1.0;
                }else{
                    v[i][j] = 0.0;
                }
                
                v[i][j] = 0;//0.1 * std::sin(pi * (x));
                       //0.1 * std::exp(- 20.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
            }
        }
    }

    template <typename T = double>
    void makeInitGradX(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& f, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& fx){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                fx[i][j] = (f[i+1][j] - f[i-1][j]) / (2.0 * dx);
            }
        }
    }

    template <typename T = double>
    void makeInitGradY(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& g, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& gy){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gy[i][j] = (g[i][j+1] - g[i][j-1]) / (2.0 * dy);
            }
        }
    }

    namespace coeff{
        template< typename T = double>
        T ax(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& ux, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind){
            if(upwind[i][j] >= 0.0){
                return (ux[i][j] + ux[i-1][j]) / (dx * dx) - 2.0 * (u[i][j] - u[i-1][j]) / (dx * dx * dx);
            }
            return (ux[i][j] + ux[i+1][j]) / (dx * dx) - 2.0 * (u[i][j] - u[i+1][j]) / (dx * dx * dx);
        }

        template< typename T = double>
        T ay(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& uy, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind){
            if(upwind[i][j] >= 0.0){
                return (uy[i][j] + uy[i][j-1]) / (dy * dy) - 2.0 * (u[i][j] - u[i][j-1]) / (dy * dy * dy);
            }
            return (uy[i][j] + uy[i][j+1]) / (dy * dy) - 2.0 * (u[i][j] - u[i][j+1]) / (dy * dy * dy);
        }
  
        template< typename T = double>
        T bx(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& ux, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind){
            if(upwind[i][j] >= 0.0){
                return 3.0 * (u[i-1][j] - u[i][j]) / (dx * dx) + (2.0 * ux[i][j] + ux[i-1][j]) / (dx);
            }
            return 3.0 * (u[i+1][j] - u[i][j]) / (dx * dx) + (2.0 * ux[i][j] + ux[i+1][j]) / (dx);
        }
        
        template< typename T = double>
        T by(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& uy, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind){
            if(upwind[i][j] >= 0.0){
                return 3.0 * (u[i][j-1] - u[i][j]) / (dy * dy) + (2.0 * uy[i][j] + uy[i][j-1]) / (dy);
            }
            return 3.0 * (u[i][j+1] - u[i][j]) / (dy * dy) + (2.0 * uy[i][j] + uy[i][j+1]) / (dy);
        }
  
        template< typename T = double>
        T cx(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& ux){
            return ux[i][j];
        }

        template< typename T = double>
        T cy(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& uy){
            return uy[i][j];
        }
            
        template< typename T = double>
        T d(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u){
            return u[i][j];
        }
    }


    template <typename T = double>
    void succX(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u1, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& ux, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(upwind[i][j] >= 0){
                    up = - upwind[i][j] * dt;
                }else{
                    up = upwind[i][j] * dt;
                }
                u2[i][j] =   coeff::ax(i, j, u1, ux, upwind) * up * up * up 
                           + coeff::bx(i, j, u1, ux, upwind) * up * up
                           + coeff::cx(i, j, ux)             * up
                           + coeff::d (i, j, u1);
            }
        }
    }

    template <typename T = double>
    void succY(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u1, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& uy, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(upwind[i][j] >= 0){
                    up = - upwind[i][j] * dt;
                }else{
                    up = upwind[i][j] * dt;
                }
                u2[i][j] =   coeff::ay(i, j, u1, uy, upwind) * up * up * up 
                           + coeff::by(i, j, u1, uy, upwind) * up * up
                           + coeff::cy(i, j, uy)             * up
                           + coeff::d (i, j, u1);
            }
        }
    }

    // ux1 -> ux2 / ux1 
    template <typename T = double>
    void succGradX(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u1, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& ux, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(upwind[i][j] >= 0){
                    up = - upwind[i][j] * dt;
                }else{
                    up = upwind[i][j] * dt;
                }
                u2[i][j] =   3.0 * coeff::ax(i, j, u1, ux, upwind) * up * up
                           + 2.0 * coeff::bx(i, j, u1, ux, upwind) * up
                           +       coeff::cx(i, j, ux);
            }
        }
    }

    
    template <typename T = double>
    void succGradY(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u1, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& uy, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& upwind, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u2){
        T up;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(upwind[i][j] >= 0){
                    up = - upwind[i][j] * dt;
                }else{
                    up = upwind[i][j] * dt;
                }
                u2[i][j] =   3.0 * coeff::ay(i, j, u1, uy, upwind) * up * up
                           + 2.0 * coeff::by(i, j, u1, uy, upwind) * up
                                 + coeff::cy(i, j, uy);
            }
        }
    }

    template <typename T = double>
    void diffusion(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u1, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u2){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u2[i][j] = u1[i][j] + dt / Re * (u1[i-1][j] - 2.0 * u1[i][j] + u1[i+1][j]) / (dx * dx) 
                                    + dt / Re * (u1[i][j-1] - 2.0 * u1[i][j] + u1[i][j+1]) / (dy * dy);
            }
        }
    }

    template <typename T = double>
    T speed(int i, int j, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& u, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& v){
        return std::sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
    }

    template <typename T = double>
    void mean(sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& f, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& g, sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>& fg){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                fg[i][j] = 0.5 * (f[i][j] + g[i][j]);
            }
        }
    }
}

int e(int i){
    int num = N + 1;
    i = ((i % num) + num) % num;
    return i;
}

int main(){

    auto u1  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto u2  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto u3  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto u4  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto ux1 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto ux2 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto ux3 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto ux4 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto uy1 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto uy2 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto uy3 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto uy4 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};

    auto v1  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto v2  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto v3  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto v4  = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vx1 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vx2 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vx3 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vx4 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vy1 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vy2 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vy3 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    auto vy4 = sksat::array_wrapper<N+1,sksat::array_wrapper<N + 1>>{};
    
    mino2357::makeInitFuncU(u1);
    mino2357::makeInitFuncV(v1);

    mino2357::makeInitGradX(u1, ux1);
    mino2357::makeInitGradY(u1, uy1);
    
    mino2357::makeInitGradX(v1, vx1);
    mino2357::makeInitGradY(v1, vy1);

    std::cout << dt * (dx * dx * Re) << std::endl;
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
            //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u1, v1));
            fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1[i][j]);
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
        
        //u
        mino2357::succX(u1, ux1, u1, u2);
        mino2357::succY(u2, uy1, v1, u3);
        //mino2357::mean(u2, u3, u3);
        mino2357::diffusion(u3, u4);

        //v
        mino2357::succX(v1, vx1, u1, v2);
        mino2357::succY(v2, vy1, v1, v3);
        //mino2357::mean(v2, v3, v3);
        mino2357::diffusion(v3, v4);

        //ux
        mino2357::succGradX(u1, ux1, u1, ux2);
        mino2357::succGradY(u1, ux2, v1, ux3);

        //uy
        mino2357::succGradX(u1, uy1, u1, uy2);
        mino2357::succGradY(u1, uy2, v1, uy3);
        
        //vx
        mino2357::succGradX(v1, vx1, u1, vx2);
        mino2357::succGradY(v1, vx2, v1, vx3);

        //vy
        mino2357::succGradX(v1, vy1, u1, vy2);
        mino2357::succGradY(v1, vy2, v1, vy3);

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                ux4[i][j] = ux3[i][j]
                            + dt / Re * (ux3[i-1][j] - 2.0 * ux3[i][j] + ux3[i+1][j]/(dx * dx))
                            + dt / Re * (uy3[i-1][j-1] + uy3[i+1][j+1] - uy3[i-1][j+1] - uy3[i+1][j-1])/(4.0 * dx * dy)
                            - dt * (ux3[i][j] * ux3[i][j] + vx3[i][j] * uy3[i][j]);
                uy4[i][j] = uy3[i][j]
                            + dt / Re * (ux3[i-1][j-1] + ux3[i+1][+1] - ux3[i-1][j+1] - ux3[i+1][j-1])/(4.0 * dx * dy)
                            + dt / Re * (uy3[i][j-1] - 2.0 * uy3[i][j] + uy3[i][j+1]/(dy * dy))
                            - dt * (uy3[i][j] * ux3[i][j] + vy3[i][j] * uy3[i][j]);
                vx4[i][j] = vx3[i][j]
                            + dt / Re * (vx3[i-1][j] - 2.0 * vx3[i][j] + vx3[i+1][j]/(dx * dx))
                            + dt / Re * (vy3[i-1][j-1] + vy3[i+1][j+1] - vy3[i-1][j+1] - vy3[i+1][j-1])/(4.0 * dx * dy)
                            - dt * (ux3[i][j] * vx3[i][j] + vx3[i][j] * vy3[i][j]);
                vy4[i][j] = vy3[i][j]
                            + dt / Re * (vx3[i-1][j-1] + vx3[i+1][j+1] - vx3[i-1][j+1] - vx3[i+1][j-1])/(4.0 * dx * dy)
                            + dt / Re * (vy3[i][j-1] - 2.0 * vy3[i][j] + vy3[i][j+1]/(dy * dy))
                            - dt * (uy3[i][j] * vx3[i][j] + vy3[i][j] * vy3[i][j]);
            }
        }

        if(it%INTV == 0){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u3, v3));
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u4[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                 u1[i][j] =  u4[i][j];
                 v1[i][j] =  v4[i][j];
                ux1[i][j] = ux4[i][j];
                uy1[i][j] = uy4[i][j];
                vx1[i][j] = vx4[i][j];
                vy1[i][j] = vy4[i][j];
            }
        }
        t += dt;
    }

    fprintf(gp, "splot '-' w l\n");
    for(int i=0; i<=N; i = i + plus){
      	for(int j=0; j<=N; j = j + plus){
           	fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1[i][j]);
      	}
       	fprintf(gp, "\n");
 	}
   	fprintf(gp, "e\n");
   	fflush(gp);

	std::cout << "simulation done.(mean)" << std::endl;
	getchar();

    pclose(gp); 
}
