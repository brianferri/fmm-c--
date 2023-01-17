#include "../lib/SMatrix.tpp"

namespace fmm
{
	/// @brief Naive Iterative Matrix multiplication : O(n^3)
	template <typename T>
	SMatrix<T> multIt(const SMatrix<T> &A, const SMatrix<T> &B);

	/// @brief Strassen Matrix multiplication : O(n^2.807)
	template <typename T>
	SMatrix<T> multSt(const SMatrix<T> &A, const SMatrix<T> &B);

	/// @brief Schonhage Laser Method Matrix multiplication : O(n^2.522) TODO
	template <typename T>
	SMatrix<T> multSl(const SMatrix<T> &A, const SMatrix<T> &B);
}

#include "../../src/fmm/fmm.ipp"