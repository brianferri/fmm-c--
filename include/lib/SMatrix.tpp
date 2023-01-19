#include <iostream>
#include <exception>
#include <initializer_list>

#define COFACTOR(i, j) minor(i, j).determinant() * ((i + j) % 2 == 0 ? 1 : -1)

/**
 * @brief Square Matrix class
 *
 * @tparam T type of data the matrix will hold
 *
 * @note the matrix is stored in a 1D array which maps to a 2D set
 * @note a matrix will have at least have O(n^2) computational complexity
 * 	for any given operation
 * 	due to the fact that it is a square matrix
 */
template <typename T>
class SMatrix
{
public:
	/* class alloc */

	SMatrix<T>(std::initializer_list<std::initializer_list<T>> list) : n(list.size()), data(new T[n * n])
	{
		if (list.size() == 0)
			throw std::invalid_argument("Matrix must have at least one row!");
		if (list.begin()->size() != list.size())
			throw std::invalid_argument("Matrix must be square!");
		size_t i = 0;
		for (auto &row : list)
		{
			size_t j = 0;
			for (auto &col : row)
			{
				data[i * n + j] = col;
				++j;
			}
			++i;
		}
	}
	SMatrix(size_t n) : n(n), data(new T[n * n]) { fill(0); }
	SMatrix() : n(0), data(nullptr) {}
	~SMatrix() { delete[] data; }

	/* util */

	size_t size() const { return n; }
	void fill(T value);
	void crop(size_t n);
	void trim(size_t n);
	SMatrix<T> *split() const;
	static SMatrix<T> merge(SMatrix<T> *m);
	void swapRows(size_t i, size_t j);
	size_t rank() const;

	/* matrix operations and algorithms */

	static SMatrix<T> identity(size_t n);
	T determinant() const;
	SMatrix<T> transpose() const;
	SMatrix<T> minor(size_t i, size_t j) const;
	SMatrix<T> adjoint() const;
	SMatrix<T> REF() const;

	/* operator overload -- should not exceed O(n^2) */
	// Implementations only consider square matrices of same size
	
	T operator()(size_t i, size_t j) const { return data[i * n + j]; }
	T &operator()(size_t i, size_t j)
	{
		if ((i >= n || j >= n) || (i < 0 || j < 0))
			throw std::range_error("Index out of range!");
		return data[i * n + j];
	}

	SMatrix<T> &operator=(const SMatrix<T> &m)
	{
		if (this == &m)
			return *this;

		delete[] data;
		n = m.size();
		data = new T[n * n];
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				data[i * n + j] = m(i, j);
			}
		}
		return *this;
	}

	friend std::ostream &operator<<(std::ostream &os, const SMatrix &m)
	{
		for (size_t i = 0; i < m.size(); ++i)
		{
			for (size_t j = 0; j < m.size(); ++j)
			{
				os << m(i, j) << " ";
			}
			os << std::endl;
		}
		return os;
	}

	friend SMatrix<T> operator+(const SMatrix<T> &A, const SMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > +");

		size_t n = A.size();
		SMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) + B(i, j);
			}
		}
		return C;
	}

	friend SMatrix<T> operator+=(SMatrix<T> &A, const SMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > +=");

		size_t n = A.size();
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				A(i, j) += B(i, j);
			}
		}
		return A;
	}

	friend SMatrix<T> operator-(const SMatrix<T> &A, const SMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > -");

		size_t n = A.size();
		SMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) - B(i, j);
			}
		}
		return C;
	}

	friend SMatrix<T> operator-=(SMatrix<T> &A, const SMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > -=");

		size_t n = A.size();
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				A(i, j) -= B(i, j);
			}
		}
		return A;
	}

	friend SMatrix<T> operator*(const SMatrix<T> &A, const T &k)
	{
		size_t n = A.size();
		SMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) * k;
			}
		}
		return C;
	}

	friend SMatrix<T> operator*(const T &k, const SMatrix<T> &A)
	{
		size_t n = A.size();
		SMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) * k;
			}
		}
		return C;
	}

protected:
	size_t n;
	T *data;
};

#include "../../src/lib/SMatrix.ipp"