#include "../lib/SquareMatrix.tpp"

namespace fmm
{
	/// @brief Naive Iterative Matrix multiplication : O(n^3)
	template <typename T>
	SquareMatrix<T> mult(const SquareMatrix<T> &A, const SquareMatrix<T> &B);

	/// @brief Strassen Matrix multiplication : O(n^2.807)
	template <typename T>
	SquareMatrix<T> strassen(const SquareMatrix<T> &A, const SquareMatrix<T> &B);

	/// @brief Schonhage Laser Method Matrix multiplication : O(n^2.522) TODO
	template <typename T>
	SquareMatrix<T> multSl(const SquareMatrix<T> &A, const SquareMatrix<T> &B);
}

#include "../../src/fmm/fmm.ipp"