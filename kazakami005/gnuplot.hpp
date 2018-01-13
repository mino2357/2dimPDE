#pragma once
#include <iostream>
#include <cmath>
#include <array>
#include "parameter.hpp"
#include "extendedArray.hpp"

using extendedArray = sksat::array_wrapper<N + 1, sksat::array_wrapper<N + 1>>;

namespace mino2357{
    class gnuplot{
    private:
        FILE* gp;
    public:
        gnuplot(){
            gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set pm3d\n");
            fprintf(gp, "set pm3d map\n");
            fprintf(gp, "set contour\n");
            //fprintf(gp, "set cbtics 0.05\n");
            fprintf(gp, "set xr [%f:%f]\n", xstart, xstart + Lx);
            fprintf(gp, "set yr [%f:%f]\n", ystart, ystart + Ly);
          	//fprintf(gp, "set zr [-0.1:1.1]\n");
            fprintf(gp, "set size square\n");
            fprintf(gp, "unset key\n");
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
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, p[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        void vectorWithLog10P(extendedArray& f, extendedArray& g, extendedArray& p){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=istart; i<=N; i = i + plus){
                for(int j=jstart; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, std::log10(std::abs(p[i][j])));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vectorWithSpeed(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, L);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
		void vectorWithSpeedLog10(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
    };
}
