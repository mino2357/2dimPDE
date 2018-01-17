/*
 * cavity flow test.
 *
 * kazakami method. SMAC method.
 * 
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>
#include <array>
#include <limits>
#include <iomanip>
#include "parameter.hpp"
#include "gnuplot.hpp"
#include "extendedArray.hpp"

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
                p[i][j] = 10.0;
            }
        }
    }

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

    template<typename T = double>
    void succDp(extendedArray& f2, extendedArray& g2, extendedArray& dp1, extendedArray& dp2, const T& t, int& itr, T& err){
        static auto R = extendedArray{};
        static auto L = extendedArray{};
        for(int k=0; ; k++){
			for(int i=1; i<N; ++i){
                for(int j=1; j<N; ++j){
                    R[i][j] = ((f2[i+1][j] - f2[i-1][j]) / (2.0 * dx) + (g2[i][j+1] - g2[i][j-1]) / (2.0 * dy)) / dt;
                    R[i][j] = R[i][j] * dx * dx * dy * dy;
                    L[i][j] = dy * dy * (dp1[i-1][j] + dp2[i+1][j]) + dx * dx * (dp1[i][j-1] + dp1[i][j+1]);
                    dp1[i][j] = - (R[i][j] - L[i][j]) / (2.0 * (dx * dx + dy * dy ));
                }
            }
            for(int i=1; i<N; ++i){
                for(int j=1; j<N; ++j){
                    dp2[i][j] = dp1[i][j];
                }
            }
			if(k%Lim == 0){
            	for(int i=1; i<N; ++i){
                	for(int j=1; j<N; ++j){
            	        R[i][j] = ((f2[i+1][j] - f2[i-1][j]) / (2.0 * dx) + (g2[i][j+1] - g2[i][j-1]) / (2.0 * dy)) / dt;
        	            R[i][j] = R[i][j] * dx * dx * dy * dy;
    	                L[i][j] = dy * dy * (dp1[i-1][j] + dp2[i+1][j]) + dx * dx * (dp1[i][j-1] + dp1[i][j+1]);
	                    err += std::abs(dp1[i][j] + (R[i][j] - L[i][j]) / (2.0 * (dx * dx + dy * dy)));
                	}
            	}
			}
            if((k > 1 && (k%Lim) == 0 && err < (Tol * N * N))){
				itr = k;
				if(t > 1.0){
					Lim = 10;
				}
				break;
            }
			err = 0.0;
        }
    }

	void setInitFuncRho(extendedArray& rho){
		double x, y;
		for(int i=0; i<=N; ++i){
			for(int j=0; j<=N; ++j){
				x = xstart + i*dx;
				y = ystart + j*dy;
				if(/*(x > 0.15 && x < 0.25) ||*/ (y > 0.45 && y < 0.55) || (x > 0.85 && x < 0.95)){// && y > 0.75 && y < 0.95){
					rho[i][j] = 1.0;
				}else{
					rho[i][j] = 0.0;
				}
			}
		}
	}
    
	void succRho(extendedArray& rho1, extendedArray& f, extendedArray& g, extendedArray& rho2){
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                if(f[i][j] >= 0.0 && g[i][j] >= 0.0){
                    rho2[i][j] = rho1[i][j] - f[i][j] * dt * (rho1[i][j] - rho1[i-1][j]) / dx
                                            - g[i][j] * dt * (rho1[i][j] - rho1[i][j-1]) / dy
                                        	+ dt * nu * (rho1[i-1][j] - 2.0 * rho1[i][j] + rho1[i+1][j]) / (dx * dx)
                                        	+ dt * nu * (rho1[i][j-1] - 2.0 * rho1[i][j] + rho1[i][j+1]) / (dy * dy);
                }else if(f[i][j] < 0.0 && g[i][j] >= 0.0){
                    rho2[i][j] = rho1[i][j] + f[i][j] * dt * (rho1[i][j] - rho1[i+1][j]) / dx
                                         	- g[i][j] * dt * (rho1[i][j] - rho1[i][j-1]) / dy
                                        	+ dt * nu * (rho1[i-1][j] - 2.0 * rho1[i][j] + rho1[i+1][j]) / (dx * dx)
                                        	+ dt * nu * (rho1[i][j-1] - 2.0 * rho1[i][j] + rho1[i][j+1]) / (dy * dy);
                }else if(f[i][j] >= 0.0 && g[i][j] < 0.0){
                    rho2[i][j] = rho1[i][j] - f[i][j] * dt * (rho1[i][j] - rho1[i-1][j]) / dx
                                        	+ g[i][j] * dt * (rho1[i][j] - rho1[i][j+1]) / dy
                                        	+ dt * nu * (rho1[i-1][j] - 2.0 * rho1[i][j] + rho1[i+1][j]) / (dx * dx)
                                        	+ dt * nu * (rho1[i][j-1] - 2.0 * rho1[i][j] + rho1[i][j+1]) / (dy * dy);
                }else{
                    rho2[i][j] = rho1[i][j] + f[i][j] * dt * (rho1[i][j] - rho1[i+1][j]) / dx
                                        	+ g[i][j] * dt * (rho1[i][j] - rho1[i][j+1]) / dy
                                        	+ dt * nu * (rho1[i-1][j] - 2.0 * rho1[i][j] + rho1[i+1][j]) / (dx * dx)
                                        	+ dt * nu * (rho1[i][j-1] - 2.0 * rho1[i][j] + rho1[i][j+1]) / (dy * dy);
                }
            }
        }
    }
}


int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 - 8);

    auto f1 = extendedArray{};
    auto f2 = extendedArray{};
    auto f3 = extendedArray{};
    
    auto g1 = extendedArray{};
    auto g2 = extendedArray{};
    auto g3 = extendedArray{};
   
    auto p1 = extendedArray{};
    auto p2 = extendedArray{};
    
	auto dp1 = extendedArray{};
    auto dp2 = extendedArray{};

	auto rho1 = extendedArray{};
	auto rho2 = extendedArray{};

	auto err = 0.0;

    mino2357::setInitFuncF(f1);
    mino2357::setInitFuncG(g1);
    mino2357::setInitFuncP(p1);
	mino2357::setInitFuncRho(rho1);

    std::cout << dt / (dx * Re) << std::endl;

    auto gp = mino2357::gnuplot();    
	gp.printRho(rho1);

    std::cout << "Enter!" << std::endl;
    getchar();

    double t = 0.0;

    for(int it=0; t<=tLimit; ++it){

		int itr = 0;

        mino2357::succF(f1, g1, p1, f2);
        mino2357::succG(g1, f1, p1, g2);
       
        mino2357::succDp(f2, g2, dp1, dp2, t, itr, err);

        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                f3[i][j] = f2[i][j] - dt * (dp2[i+1][j] - dp2[i-1][j]) / (2.0 * dx);
                g3[i][j] = g2[i][j] - dt * (dp2[i][j+1] - dp2[i][j-1]) / (2.0 * dy);
            }
        }

        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                p2[i][j] = p1[i][j] + dp2[i][j];
            }
        }

		mino2357::succRho(rho1, f3, g3, rho2);

        if(it%INTV == 0){
            std::cout << itr << "  " << t << "  " << err << std::endl;
            //gp.printP(p2);
            //gp.speed(f3, g3);
            //gp.vector(f3, g3);
            //gp.vectorWithP(f3, g3, p2);
            //gp.vectorWithSpeed(f3, g3);
            //gp.vectorWithSpeedLog10(f3, g3);
            //gp.vectorWithLog10P(f3, g3, p2);
            //gp.multiplot(f3, g3, p2, t);
            //gp.multiplotMakePNG(f3, g3, p2, t);
			//gp.printRho(rho2);	
			gp.multiplotMakePNGWithDensity(f3, g3, rho2, t);
        }
   
        for(int i=1; i<N; ++i){
            for(int j=1; j<N; ++j){
                f1[i][j] = f3[i][j];
                g1[i][j] = g3[i][j];
                p1[i][j] = p2[i][j];
				rho1[i][j] = rho2[i][j];
            }
        }
		//境界圧力
        for(int ij=0; ij<=N; ++ij){
            p1[ij][0] = p1[ij][1];
            p1[0][ij] = p1[1][ij];
            p1[N][ij] = p1[N][ij-1];
        }
		//cavity flow
        for(int i=0; i<=N; ++i){
            f1[i][N] = 1.0;
        }
        t += dt;
		p1[1][1] = 10.0;
    }
    gp.print(f1);
}
