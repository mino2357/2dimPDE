/*
 * guruguru test.
 *
 * kazakami method.
 * 
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>
#include <array>

constexpr int N = 128;
constexpr double Re = 1000;
constexpr double Lx =  1.0;
constexpr double Ly =  1.0;
constexpr double xstart = 0.0;
constexpr double ystart = 0.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.005;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 2000000.0 * pi;
constexpr double Tol = 10e-10; 
constexpr int INTV = 5;
constexpr int plus = 2;

namespace rittai3d{
	namespace utility{
		// [ minimum, maximum ) の範囲でラップアラウンド
		template <typename T>
		constexpr T wrap_around(T value, T minimum, T maximum){
			const T n = (value - minimum) % (maximum - minimum);
			return n >= 0 ? (n + minimum) : (n + maximum); 
		}
	}
}

namespace sksat {
	template<std::size_t Num, typename T = double>
	class array_wrapper{
	private:
		std::array<T, Num> arr;
	public:
		constexpr array_wrapper() : arr() {}
		~array_wrapper() = default;

		T& operator[](int i){
			i = rittai3d::utility::wrap_around(i, 0, static_cast<int>(Num));
			return arr[i];
		}
	};
}

using extendedArray = sksat::array_wrapper<N + 1, sksat::array_wrapper<N + 1>>;

namespace mino2357{
    
    template <typename T = double>
    void setInitFuncF(extendedArray& f){
        T x, y;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.2 && x <= 0.6 && y >= 0.2 && y <= 0.6){
                    f[i][j] = 0.0;
                }else{
                    f[i][j] = 0.0;
                }
            } 
        }
        for(int i=0; i<=N; ++i){
            f[i][0] = 0.0;
            f[i][N] = 1.0;
        }
    }
    
    template <typename T = double>
    void setInitFuncG(extendedArray& g){
        T x, y;
        //T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.2 && x <= 0.6 && y >= 0.2 && y <= 0.6){
                    g[i][j] = 0.0;
                }else{
                    g[i][j] = 0.0;
                }
            } 
        }
        for(int j=0; j<=N; ++j){
            g[0][j] = 0.0;
            g[N][j] = 0.0;
        }
    }

    template<typename T = double>
    void setInitFuncP(extendedArray& p){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                p[i][j] = 0.0;
            }
        }
    }
    

    class gnuplot{
    private:
        FILE* gp;
    public:
        gnuplot(){
            gp = popen("gnuplot -persist", "w");
            //fprintf(gp, "set pm3d\n");
            //fprintf(gp, "set pm3d map\n");
            //fprintf(gp, "set contour\n");
            //fprintf(gp, "set cbtics 0.05\n");
            fprintf(gp, "set xr [%f:%f]\n", xstart, xstart + Lx);
            fprintf(gp, "set yr [%f:%f]\n", ystart, ystart + Ly);
          	//fprintf(gp, "set zr [-0.1:1.1]\n");
            fprintf(gp, "set size square\n");
        }
        ~gnuplot(){ pclose(gp);}

        void print(extendedArray& f){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u1, v1));
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, f[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void speed(extendedArray& f, extendedArray& g){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u1, v1));
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void plot(double x, double y){
            fprintf(gp, "plot '-' w l\n");
            fprintf(gp, "%f %f\n", x, y);
            fprintf(gp, "\n");
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vector(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j], g[i][j],
                                                    std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vectorWithP(extendedArray& f, extendedArray& g, extendedArray& p){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = 20 * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, p[i][j]);
                    //fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j], g[i][j], p[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
    };

    template <typename T = double>
    void succF(extendedArray& f1, extendedArray& g1, extendedArray& p, extendedArray& f2){
        //for(int i=0; i<=N; ++i){
        //    for(int j=0; j<=N; ++j){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                if(f1[i][j] >= 0.0 && g1[i][j] >= 0.0){
                    f2[i][j] = f1[i][j] - dt * (p[i+1][j] - p[i-1][j]) / (2.0 * dx)
                                        - f1[i][j] * dt * (f1[i][j] - f1[i-1][j]) / dx
                                        - g1[i][j] * dt * (f1[i][j] - f1[i][j-1]) / dy
                                        + dt / Re * (f1[i-1][j] - 2.0 * f1[i][j] + f1[i+1][j]) / (dx * dx)
                                        + dt / Re * (f1[i][j-1] - 2.0 * f1[i][j] + f1[i][j+1]) / (dy * dy);
                }else if(f1[i][j] < 0.0 && g1[i][j] >= 0.0){
                    f2[i][j] = f1[i][j] - dt * (p[i+1][j] - p[i-1][j]) / (2.0 * dx)
                                        + f1[i][j] * dt * (f1[i][j] - f1[i+1][j]) / dx
                                        - g1[i][j] * dt * (f1[i][j] - f1[i][j-1]) / dy
                                        + dt / Re * (f1[i-1][j] - 2.0 * f1[i][j] + f1[i+1][j]) / (dx * dx)
                                        + dt / Re * (f1[i][j-1] - 2.0 * f1[i][j] + f1[i][j+1]) / (dy * dy);
                }else if(f1[i][j] >= 0.0 && g1[i][j] < 0.0){
                    f2[i][j] = f1[i][j] - dt * (p[i+1][j] - p[i-1][j]) / (2.0 * dx)
                                        - f1[i][j] * dt * (f1[i][j] - f1[i-1][j]) / dx
                                        + g1[i][j] * dt * (f1[i][j] - f1[i][j+1]) / dy
                                        + dt / Re * (f1[i-1][j] - 2.0 * f1[i][j] + f1[i+1][j]) / (dx * dx)
                                        + dt / Re * (f1[i][j-1] - 2.0 * f1[i][j] + f1[i][j+1]) / (dy * dy);
                }else{
                    f2[i][j] = f1[i][j] - dt * (p[i+1][j] - p[i-1][j]) / (2.0 * dx)
                                        + f1[i][j] * dt * (f1[i][j] - f1[i+1][j]) / dx
                                        + g1[i][j] * dt * (f1[i][j] - f1[i][j+1]) / dy
                                        + dt / Re * (f1[i-1][j] - 2.0 * f1[i][j] + f1[i+1][j]) / (dx * dx)
                                        + dt / Re * (f1[i][j-1] - 2.0 * f1[i][j] + f1[i][j+1]) / (dy * dy);
                }
            }
        }
    }
    
    template <typename T = double>
    void succG(extendedArray& g1, extendedArray& f1, extendedArray& p, extendedArray& g2){
        //for(int i=0; i<=N; ++i){
        //    for(int j=0; j<=N; ++j){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                if(f1[i][j] >= 0.0 && g1[i][j] >= 0.0){
                    g2[i][j] = g1[i][j] - dt * (p[i][j+1] - p[i][j-1]) / (2.0 * dy)
                                        - f1[i][j] * dt * (g1[i][j] - g1[i-1][j]) / dx
                                        - g1[i][j] * dt * (g1[i][j] - g1[i][j-1]) / dy
                                        + dt / Re * (g1[i-1][j] - 2.0 * g1[i][j] + g1[i+1][j]) / (dx * dx)
                                        + dt / Re * (g1[i][j-1] - 2.0 * g1[i][j] + g1[i][j+1]) / (dy * dy);
                }else if(f1[i][j] < 0.0 && g1[i][j] >= 0.0){
                    g2[i][j] = g1[i][j] - dt * (p[i][j+1] - p[i][j-1]) / (2.0 * dy)
                                        + f1[i][j] * dt * (g1[i][j] - g1[i+1][j]) / dx
                                        - g1[i][j] * dt * (g1[i][j] - g1[i][j-1]) / dy
                                        + dt / Re * (g1[i-1][j] - 2.0 * g1[i][j] + g1[i+1][j]) / (dx * dx)
                                        + dt / Re * (g1[i][j-1] - 2.0 * g1[i][j] + g1[i][j+1]) / (dy * dy);
                }else if(f1[i][j] >= 0.0 && g1[i][j] < 0.0){
                    g2[i][j] = g1[i][j] - dt * (p[i][j+1] - p[i][j-1]) / (2.0 * dy)
                                        - f1[i][j] * dt * (g1[i][j] - g1[i-1][j]) / dx
                                        + g1[i][j] * dt * (g1[i][j] - g1[i][j+1]) / dy
                                        + dt / Re * (g1[i-1][j] - 2.0 * g1[i][j] + g1[i+1][j]) / (dx * dx)
                                        + dt / Re * (g1[i][j-1] - 2.0 * g1[i][j] + g1[i][j+1]) / (dy * dy);
                }else{
                    g2[i][j] = g1[i][j] - dt * (p[i][j+1] - p[i][j-1]) / (2.0 * dy)
                                        + f1[i][j] * dt * (g1[i][j] - g1[i+1][j]) / dx
                                        + g1[i][j] * dt * (g1[i][j] - g1[i][j+1]) / dy
                                        + dt / Re * (g1[i-1][j] - 2.0 * g1[i][j] + g1[i+1][j]) / (dx * dx)
                                        + dt / Re * (g1[i][j-1] - 2.0 * g1[i][j] + g1[i][j+1]) / (dy * dy);
                }
            }
        }
    }
}


int main(){

    auto f1 = extendedArray{};
    auto f2 = extendedArray{};
    auto f3 = extendedArray{};
    
    auto g1 = extendedArray{};
    auto g2 = extendedArray{};
    auto g3 = extendedArray{};
   
    auto p1 = extendedArray{};
    auto dp2 = extendedArray{};
    auto dp3 = extendedArray{};
    auto p4 = extendedArray{};
    
    auto R = extendedArray{};
    auto L = extendedArray{};

    mino2357::setInitFuncF(f1);
    mino2357::setInitFuncG(g1);
    mino2357::setInitFuncP(p1);

    std::cout << dt / (dx * Re) << std::endl;

    auto gp = mino2357::gnuplot();    

    gp.vectorWithP(f1, g2, p1);
    //gp.vector(f1, g2);
    //gp.print(p1);
    std::cout << "Enter!" << std::endl;
    getchar();

    double t = 0.0;

    for(int it=0; t<=tLimit; ++it){

        mino2357::succF(f1, g1, p1, f2);
        mino2357::succG(g1, f1, p1, g2);
       
        //int times = 10000;
 
        double a,b;

        for(int k=0; ; k++){
            //for(int i=0; i<=N; ++i){
            //    for(int j=0; j<=N; ++j){
            //if(k == times-2){
                a = dp2[N-10][10];
            //}
            for(int i=1; i<N; ++i){
                for(int j=1; j<N; ++j){
                    R[i][j] = ((f2[i+1][j] - f2[i-1][j]) / (2.0 * dx) + (g2[i][j+1] - g2[i][j-1]) / (2.0 * dy)) / dt;
                    R[i][j] = R[i][j] * dx * dx * dy * dy;
                    L[i][j] = dy * dy * (dp2[i-1][j] + dp2[i+1][j])  + dx * dx * (dp2[i][j-1] + dp2[i][j+1]);
                    dp3[i][j] = - (R[i][j] - L[i][j]) / (2.0 * (dx * dx + dy * dy ));
                }
            }
            for(int j=0; j<=N; ++j){
                dp3[0][j] = dp3[1][j];
                dp3[N][j] = dp3[N-1][j];
            }
            for(int i=0; i<=N; ++i){
                dp3[i][0] = dp3[i][1];
            }
            //if(k == times-1){
                b = dp3[N-10][10];
            //}
            if(std::abs((a-b)/a) < Tol){
                break;
            }
            //for(int i=0; i<=N; ++i){
            //    for(int j=0; j<=N; ++j){
            for(int i=1; i<N; ++i){
                for(int j=1; j<N; ++j){
                    dp2[i][j] = dp3[i][j];
                }
            }
        }

        //for(int i=0; i<=N; ++i){
        //    for(int j=0; j<=N; ++j){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                f3[i][j] = f2[i][j] - dt * (dp3[i+1][j] - dp3[i-1][j]) / (2.0 * dx);
                g3[i][j] = g2[i][j] - dt * (dp3[i][j+1] - dp3[i][j-1]) / (2.0 * dy);
            }
        }

        //for(int i=0; i<=N; ++i){
        //    for(int j=0; j<=N; ++j){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                p4[i][j] = p1[i][j] + dp3[i][j];
            }
        }

        if(it%INTV == 0){
            std::cout << t << " " << std::abs((a-b)/a) << std::endl;
            //gp.print(f3);
            //gp.print(g3);
            //gp.print(p4);
            //gp.speed(f3, g3);
            //gp.vector(f3, g3);
            gp.vectorWithP(f3, g3, p4);
        }
   
        //for(int i=0; i<=N; ++i){
        //    for(int j=0; j<=N; ++j){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                f1[i][j] = f3[i][j];
                g1[i][j] = g3[i][j];
                p1[i][j] = p4[i][j];
                dp2[i][j] = 0.0;
            }
        }
        for(int i=0; i<=N; ++i){
            f1[i][N] = 1.0;
        }
        t += dt;
    }
    gp.print(f1);
}
