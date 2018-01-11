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

constexpr int N = 256;
constexpr double Lx =  2.0;
constexpr double Ly =  2.0;
constexpr double xstart = - 1.0;
constexpr double ystart = - 1.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = 0.001;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 2.0 * pi;
constexpr int INTV = 200;
constexpr int plus = 4;

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
    void makeInitFuncF(extendedArray& f){
        T x, y;
        //T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.2 && x <= 0.6 && y >= 0.2 && y <= 0.6){
                    f[i][j] = 1.0;
                }else{
                    f[i][j] = 0.0;
                }
                
               f[i][j] = (std::sin(pi * (x + y)));// + std::sin(pi * (y));
                       //0.1 * std::exp(- 20.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
                
            } 
        }
    }
    
    template <typename T = double>
    void makeInitFuncG(extendedArray& g){
        T x, y;
        //T a = 0.0; T b = 0.0;
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                x = xstart + i * dx;
                y = ystart + j * dy;
                
                if(x >= 0.2 && x <= 0.6 && y >= 0.2 && y <= 0.6){
                    g[i][j] = 1.0;
                }else{
                    g[i][j] = 0.0;
                }
                
               //f[i][j] = //.1 * (std::sin(pi * (x)));// + std::sin(pi * (y));
                       //0.1 * std::exp(- 20.0 * (std::pow(x - a, 2) + std::pow(y - b, 2)));
                
            } 
        }
    }

    

    class gnuplot{
    private:
        FILE* gp;
    public:
        gnuplot(){
            gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set pm3d\n");
            //fprintf(gp, "set pm3d map\n");
            fprintf(gp, "set contour\n");
            fprintf(gp, "set xr [%f:%f]\n", xstart, xstart + Lx);
            fprintf(gp, "set yr [%f:%f]\n", ystart, ystart + Ly);
          	//fprintf(gp, "set zr [-0.1:1.1]\n");
            fprintf(gp, "set size square\n");
        }
        ~gnuplot(){ pclose(gp);}

        void print(extendedArray f){
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
    };

    template <typename T = double>
    void succ(extendedArray& f1, extendedArray& g1,extendedArray& f2){
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                if(f1[i][j] >= 0.0 && g1[i][j] >= 0.0){
                    f2[i][j] = f1[i][j] - f1[i][j] * dt * (f1[i][j] - f1[i-1][j]) / dx
                                        - g1[i][j] * dt * (f1[i][j] - f1[i][j-1]) / dy;
                }else if(f1[i][j] < 0.0 && g1[i][j] >= 0.0){
                    f2[i][j] = f1[i][j] + f1[i][j] * dt * (f1[i][j] - f1[i+1][j]) / dx
                                        - g1[i][j] * dt * (f1[i][j] - f1[i][j-1]) / dy;
                }else if(f1[i][j] >= 0.0 && g1[i][j] < 0.0){
                    f2[i][j] = f1[i][j] - f1[i][j] * dt * (f1[i][j] - f1[i-1][j]) / dx
                                        + g1[i][j] * dt * (f1[i][j] - f1[i][j+1]) / dy;
                }else{
                    f2[i][j] = f1[i][j] + f1[i][j] * dt * (f1[i][j] - f1[i+1][j]) / dx
                                        + g1[i][j] * dt * (f1[i][j] - f1[i][j+1]) / dy;
                }
            }
        }
    }

}


int main(){

    auto f1 = extendedArray{};
    auto f2 = extendedArray{};
    
    auto g1 = extendedArray{};
    auto g2 = extendedArray{};
    
    mino2357::makeInitFuncF(f1);
    mino2357::makeInitFuncG(g1);

    std::cout << dt / dx << std::endl;

    auto gp = mino2357::gnuplot();    

    gp.print(f1);
    std::cout << "Enter!" << std::endl;
    getchar();


    double t = 0.0;

    for(int it=0; t<=tLimit; ++it){

        mino2357::succ(f1, g1, f2);
        mino2357::succ(g1, f1, g2);

        if(it%INTV == 0){
            gp.print(f2);
        }
   
        for(int i=0; i<=N; ++i){
            for(int j=0; j<=N; ++j){
                f1[i][j] = f2[i][j];
                g1[i][j] = g2[i][j];
            }
        }
        t += dt;
    }
    gp.print(f1);
}
