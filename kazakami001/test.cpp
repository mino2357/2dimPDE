#include <iostream>
#include <array>

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
}

namespace sksat {
	template<std::size_t N, typename T = double>
	class array_wrapper {
	private:
		std::array<T,N> arr;
	public:
		constexpr array_wrapper() : arr() {}
		~array_wrapper() = default;

		T& operator[](int i){
			i = mino2357::utility::wrap_around(i,0,static_cast<int>(N));
			return arr[i];
		}
	};
}

int main(int argc, char **argv){
	constexpr int N = 2;
/*	sksat::array_wrapper<N> test;
	test[0] = 1.0;
	test[1] = 2.0;
	test[2] = 3.0;
	for(auto i=0;i<(N*2);i++){
		std::cout<<i<<": "<<test[i]<<std::endl;
	}
*/
	sksat::array_wrapper<N,sksat::array_wrapper<N>> a;
	a[0][0] = 1.0;
	a[0][1] = 2.0;
	a[1][0] = 3.0;
	a[1][1] = 4.0;
	for(int i=-1;i<=N+1;i++){
		for(auto j=-1;j<=N+1;j++)
			std::cout<<"a["<<i<<"]["<<j<<"]: "<<a[i][j]<<std::endl;
	}
	return 0;
}
