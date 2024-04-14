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
 * 	for any given operation (except lookup)
 * 	due to the fact that it is a square matrix
 */
template <typename T>
class SquareMatrix
{
public:
	/* class alloc */

	SquareMatrix(size_t n) : n(n), data(new T[n * n]) {}
	SquareMatrix(size_t n, T value) : n(n), data(new T[n * n]) { fill(value); }
	SquareMatrix() : n(0), data(nullptr) {}
	SquareMatrix(std::initializer_list<std::initializer_list<T>> list) : n(list.size()), data(new T[n * n])
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
	SquareMatrix(const SquareMatrix<T> &m) : n(m.size()), data(new T[n * n])
	{
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				data[i * n + j] = m(i, j);
			}
		}
	}
	~SquareMatrix() { delete[] data; }

	/* util */

	size_t size() const { return n; }
	void fill(T value);
	void crop(size_t n);
	void trim(size_t n);
	SquareMatrix<T> *split() const;
	static SquareMatrix<T> merge(SquareMatrix<T> *m);
	void swapRows(size_t i, size_t j);
	size_t rank() const;

	/* matrix operations and algorithms */

	static SquareMatrix<T> identity(size_t n);
	static SquareMatrix<T> kronecker(const SquareMatrix<T> &A, const SquareMatrix<T> &B);
	T determinant() const;
	SquareMatrix<T> transpose() const;
	SquareMatrix<T> minor(size_t i, size_t j) const;
	SquareMatrix<T> adjoint() const;
	SquareMatrix<T> REF() const;
	SquareMatrix<T> RREF() const;
	SquareMatrix<T> inverse() const;

	/* operator overload -- should not exceed O(n^2) */
	// Implementations only consider square matrices of same size
	T operator()(size_t i, size_t j) const { return data[i * n + j]; }
	T &operator()(size_t i, size_t j)
	{
		if (i >= n || j >= n)
			throw std::range_error("Index out of range!");
		return data[i * n + j];
	}
	T operator[](size_t i) const { return data[i]; }
	T &operator[](size_t i)
	{
		if (i >= n * n)
			throw std::range_error("Index out of range!");
		return data[i];
	}

	SquareMatrix<T> &operator=(const SquareMatrix<T> &m)
	{
		if (this == &m) return *this;
		if (n != m.size())
		{
			delete[] data;
			n = m.size();
			data = new T[n * n];
		}
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				data[i * n + j] = m(i, j);
			}
		}
		return *this;
	}

	friend std::ostream &operator<<(std::ostream &os, const SquareMatrix &m)
	{
		os << "[";
		for (size_t i = 0; i < m.size(); ++i)
		{
			if (i > 0)
				os << " ";
			os << "[";
			for (size_t j = 0; j < m.size(); ++j)
			{
				os << m(i, j);
				if (j < m.size() - 1)
					os << ", ";
			}
			os << "]";
			if (i < m.size() - 1)
				os << std::endl;
		}
		os << "]";
		return os;
	}

	friend SquareMatrix<T> operator+(const SquareMatrix<T> &A, const SquareMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > +");

		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) + B(i, j);
			}
		}
		return C;
	}

	friend SquareMatrix<T> operator+=(SquareMatrix<T> &A, const SquareMatrix<T> &B)
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

	friend SquareMatrix<T> operator-(const SquareMatrix<T> &A, const SquareMatrix<T> &B)
	{
		if (A.size() != B.size())
			throw std::range_error("Sizes are incompatible! > -");

		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) - B(i, j);
			}
		}
		return C;
	}

	friend SquareMatrix<T> operator-=(SquareMatrix<T> &A, const SquareMatrix<T> &B)
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

	friend SquareMatrix<T> operator*(const SquareMatrix<T> &A, const T &k)
	{
		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) * k;
			}
		}
		return C;
	}

	friend SquareMatrix<T> operator*(const T &k, const SquareMatrix<T> &A)
	{
		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) * k;
			}
		}
		return C;
	}

	friend SquareMatrix<T> operator/(const SquareMatrix<T> &A, const T &k)
	{
		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) / k;
			}
		}
		return C;
	}

	friend SquareMatrix<T> operator/(const T &k, const SquareMatrix<T> &A)
	{
		size_t n = A.size();
		SquareMatrix<T> C(n);
		for (size_t i = 0; i < n; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) = A(i, j) / k;
			}
		}
		return C;
	}

protected:
	size_t n;
	T *data;
};

#include "../../src/lib/SquareMatrix.ipp"