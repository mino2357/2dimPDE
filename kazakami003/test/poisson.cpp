/*
 * poisson test.
 *
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>
#include <array>
#include <Eigen/Sparse>

constexpr int N = 10;
/*
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
*/
int main(){

    Eigen::SparseMatrix<double> A(N, N);
    Eigen::VectorXd b(N), x(N);
    
    for(int i=0; i<N; ++i){
        b(i) = 0.0;
    }

    b(0)   = 0.0;
    b(N-1) = - 1.0;

    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            if(i == j){
                A.insert(i,j) = - 2.0;
            }else if(i == (j+1) || i == (j-1)){
                A.insert(i,j) = 1.0;
            }else{
                A.insert(i,j) = 0.0;
            }
        }
    }
    
    std::cout << A << std::endl;
    std::cout << std::endl;
    std::cout << b << std::endl;
    std::cout << std::endl;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(A);
    x = cg.solve(b);

    std::cout << x << std::endl;

}
