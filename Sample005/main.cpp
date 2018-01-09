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

constexpr int N = 64;
constexpr double Re = 100.0;
constexpr double Lx =  2.0;
constexpr double Ly =  2.0;
constexpr double xstart = -1.0;
constexpr double ystart = -1.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.0001;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 20000.0 * pi;
constexpr int INTV = 100;
constexpr int plus = 2;

namespace mino2357{
    namespace utility{
		// [ minimum, maximum ) の範囲でラップアラウンド
		template <typename T>
		constexpr T wrap_around( T value, T minimum, T maximum )
		{
			const T n = ( value - minimum ) % ( maximum - minimum );
			return n >= 0 ? ( n + minimum ) : ( n + maximum ); 
		}
	}
	
    template <std::size_t N, typename T = double>
    class extendedArray{
	
		std::array< std::array< T, N >, N > arr_;
	
    public:
        constexpr extendedArray() : arr_() {
        }
		
        ~extendedArray() = default;

		auto& operator[]( int i ) {
			i = utility::wrap_around( i, 0, static_cast< int >( N ) );
			return arr_[i];
		}
		
		constexpr T& operator()(int i, int j) {
            i = utility::wrap_around( i, 0, static_cast< int >( N ) );
			j = utility::wrap_around( j, 0, static_cast< int >( N ) );
            return arr_[ i ][ j ];
		}
    };
    
    template <typename T = double>
    void makeInitFuncU(extendedArray<N + 1>& u){
        T x, y;
        T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.2 && x <= 0.6 && y >= 0.2 && y <= 0.6){
                    u[i][j] = 1.0;
                }else{
                    u[i][j] = 0.0;
                }
                /*
                u[i][j] = //0.2 * std::sin(pi * (x));
                          2.0 * std::exp(- 10.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
                */
            }
        }
    }

    template <typename T = double>
    void makeInitFuncV(extendedArray<N + 1>& v){
        
        T x, y;
        T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
              /*
                if(xstart+i*dx >= 0.2 && xstart+i*dx <= 0.6 && ystart+j*dy >= 0.2 && ystart+j*dy <= 0.6){
                    v[i][j] = 1.0;
                }else{
                    v[i][j] = 0.0;
                }
                */
                v[i][j] = 0.0;//std::sin(pi * (x + y)) + 2.0;
                         //- 2.0 * std::exp(- 10.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
            }
        }
    }

    template <typename T = double>
    void makeInitGradX(extendedArray<N + 1>& f, extendedArray<N + 1>& fx){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                fx[i][j] = (f[i+1][j] - f[i-1][j]) / (2.0 * dx);
            }
        }
    }

    template <typename T = double>
    void makeInitGradY(extendedArray<N + 1>& g, extendedArray<N + 1>& gy){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                gy[i][j] = (g[i][j+1] - g[i][j-1]) / (2.0 * dy);
            }
        }
    }

    namespace coeff{
        template< typename T = double>
        T ax(int i, int j, extendedArray<N + 1>& u, extendedArray<N + 1>& ux, extendedArray<N + 1>& upwind){
            if(upwind[i][j] >= 0.0){
                return (ux[i][j] + ux[i-1][j]) / (dx * dx) - 2.0 * (u[i][j] - u[i-1][j]) / (dx * dx * dx);
            }
            return (ux[i][j] + ux[i+1][j]) / (dx * dx) - 2.0 * (u[i][j] - u[i+1][j]) / (dx * dx * dx);
        }

        template< typename T = double>
        T ay(int i, int j, extendedArray<N + 1>& u, extendedArray<N + 1>& uy, extendedArray<N + 1>& upwind){
            if(upwind[i][j] >= 0.0){
                return (uy[i][j] + uy[i][j-1]) / (dx * dx) - 2.0 * (u[i][j] - u[i][j-1]) / (dx * dx * dx);
            }
            return (uy[i][j] + uy[i][j+1]) / (dx * dx) - 2.0 * (u[i][j] - u[i][j+1]) / (dx * dx * dx);
        }
  
        template< typename T = double>
        T bx(int i, int j, extendedArray<N + 1>& u, extendedArray<N + 1>& ux, extendedArray<N + 1>& upwind){
            if(upwind[i][j] >= 0.0){
                return 3.0 * (u[i-1][j] - u[i][j]) / (dx * dx) + (2.0 * ux[i][j] + ux[i-1][j]) / (dx);
            }
            return 3.0 * (u[i+1][j] - u[i][j]) / (dx * dx) + (2.0 * ux[i][j] + ux[i+1][j]) / (dx);
        }
        
        template< typename T = double>
        T by(int i, int j, extendedArray<N + 1>& u, extendedArray<N + 1>& uy, extendedArray<N + 1>& upwind){
            if(upwind[i][j] >= 0.0){
                return 3.0 * (u[i][j-1] - u[i][j]) / (dx * dx) + (2.0 * uy[i][j] + uy[i][j-1]) / (dx);
            }
            return 3.0 * (u[i][j+1] - u[i][j]) / (dx * dx) + (2.0 * uy[i][j] + uy[i][j+1]) / (dx);
        }
  
        template< typename T = double>
        T cx(int i, int j, extendedArray<N + 1>& ux){
            return ux[i][j];
        }

        template< typename T = double>
        T cy(int i, int j, extendedArray<N + 1>& uy){
            return uy[i][j];
        }
            
        template< typename T = double>
        T d(int i, int j, extendedArray<N + 1>& u){
            return u[i][j];
        }
    }


    template <typename T = double>
    void succX(extendedArray<N + 1>& u1, extendedArray<N + 1>& ux, extendedArray<N + 1>& upwind, extendedArray<N + 1>& u2){
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

    // ux1 -> ux2 / ux1 
    template <typename T = double>
    void succGradX(extendedArray<N + 1>& u1, extendedArray<N + 1>& ux, extendedArray<N + 1>& upwind, extendedArray<N + 1>& u2){
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
    void succY(extendedArray<N + 1>& u1, extendedArray<N + 1>& uy, extendedArray<N + 1>& upwind, extendedArray<N + 1>& u2){
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

    template <typename T = double>
    void succGradY(extendedArray<N + 1>& u1, extendedArray<N + 1>& uy, extendedArray<N + 1>& upwind, extendedArray<N + 1>& u2){
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
    void diffusion(extendedArray<N + 1>& u1, extendedArray<N + 1>& u2){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                u2[i][j] = u2[i][j] + dt / Re * (u1[i-1][j] - 2.0 * u1[i][j] + u1[i+1][j]) / (dx * dx) 
                                    + dt / Re * (u1[i][j-1] - 2.0 * u1[i][j] + u1[i][j+1]) / (dy * dy);
            }
        }
    }

    template <typename T = double>
    T speed(int i, int j, extendedArray<N + 1>& u, extendedArray<N + 1>& v){
        return std::sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
    }

    template <typename T = double>
    void mean(extendedArray<N + 1>& f, extendedArray<N + 1> g, extendedArray<N + 1> fg){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                fg[i][j] = 0.5 * (f[i][j] + g[i][j]);
            }
        }
    }
}


int main(){

    auto u1  = mino2357::extendedArray<N + 1>{};
    auto u2  = mino2357::extendedArray<N + 1>{};
    auto u3  = mino2357::extendedArray<N + 1>{};
    auto ux1 = mino2357::extendedArray<N + 1>{};
    auto ux2 = mino2357::extendedArray<N + 1>{};
    auto ux3 = mino2357::extendedArray<N + 1>{};
    auto uy1 = mino2357::extendedArray<N + 1>{};
    auto uy2 = mino2357::extendedArray<N + 1>{};
    auto uy3 = mino2357::extendedArray<N + 1>{};

    auto v1  = mino2357::extendedArray<N + 1>{};
    auto v2  = mino2357::extendedArray<N + 1>{};
    auto v3  = mino2357::extendedArray<N + 1>{};
    auto vx1 = mino2357::extendedArray<N + 1>{};
    auto vx2 = mino2357::extendedArray<N + 1>{};
    auto vx3 = mino2357::extendedArray<N + 1>{};
    auto vy1 = mino2357::extendedArray<N + 1>{};
    auto vy2 = mino2357::extendedArray<N + 1>{};
    auto vy3 = mino2357::extendedArray<N + 1>{};

    mino2357::makeInitFuncU(u1);
    mino2357::makeInitFuncV(v1);

    mino2357::makeInitGradX(u1, ux1);
    mino2357::makeInitGradY(u1, uy1);
    
    mino2357::makeInitGradX(v1, vx1);
    mino2357::makeInitGradY(v1, vy1);

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
            fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u1[i][j]);//mino2357::speed(i, j, u1, v1));
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
        mino2357::diffusion(u1, u3);

        //v
        mino2357::succX(v1, vx1, u1, v2);
        mino2357::succY(v2, vy1, v1, v2);
        //mino2357::mean(v2, v3, v3);
        mino2357::diffusion(v1, v3);

        //ux
        mino2357::succGradX(ux1, ux1, u1, ux2);
        mino2357::succGradY(ux2, ux1, v1, ux3);

        //uy
        mino2357::succGradX(vx1, uy1, u1, vx2);
        mino2357::succGradY(vx2, uy1, v1, vx3);
        
        //vx
        mino2357::succGradX(uy1, vx1, u1, uy2);
        mino2357::succGradY(uy2, vx1, v1, uy3);

        //vy
        mino2357::succGradX(vy1, vy1, u1, vy2);
        mino2357::succGradY(vy2, vy1, v1, vy3);

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                ux3[i][j] = ux3[i][j] - dt * (ux1[i][j] * ux1[i][j] + vx1[i][j] * uy1[i][j]);
                uy3[i][j] = uy3[i][j] - dt * (uy1[i][j] * ux1[i][j] + vy1[i][j] * uy1[i][j]);
                vx3[i][j] = vx3[i][j] - dt * (ux1[i][j] * vx1[i][j] + vx1[i][j] * vy1[i][j]);
                vy3[i][j] = vy3[i][j] - dt * (uy1[i][j] * vx1[i][j] + vy1[i][j] * vy1[i][j]);
            }
        }

        if(it%INTV == 0){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, u3(i, j));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                 u1[i][j] =  u3[i][j];
                 v1[i][j] =  v3[i][j];
                ux1[i][j] = ux3[i][j];
                uy1[i][j] = uy3[i][j];
                vx1[i][j] = vx3[i][j];
                vy1[i][j] = vy3[i][j];
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
