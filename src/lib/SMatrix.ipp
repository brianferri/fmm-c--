template <typename T>
void SMatrix<T>::fill(T value)
{
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			data[i * n + j] = value;
		}
	}
}

template <typename T>
void SMatrix<T>::crop(size_t n)
{
	if (n < this->n)
		throw std::invalid_argument("If you want to trim, use trim()!");
	if (n == this->n)
		return;

	SMatrix<T> tmp(n);
	tmp.fill(0);
	for (size_t i = 0; i < this->n; ++i)
	{
		for (size_t j = 0; j < this->n; ++j)
		{
			tmp(i, j) = data[i * this->n + j];
		}
	}
	delete[] data;
	data = tmp.data;
	tmp.data = nullptr;
	this->n = n;
}

template <typename T>
void SMatrix<T>::trim(size_t n)
{
	if (n > this->n)
		throw std::invalid_argument("If you want to crop, use crop()!");
	if (n == this->n)
		return;

	SMatrix<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			tmp(i, j) = data[i * this->n + j];
		}
	}
	delete[] data;
	data = tmp.data;
	tmp.data = nullptr;
	this->n = n;
}

template <typename T>
SMatrix<T> *SMatrix<T>::split() const
{
	size_t n = this->n % 2 == 0 ? this->n / 2 : this->n / 2 + 1;
	SMatrix<T> *m = new SMatrix<T>[4];
	for (size_t i = 0; i < 4; ++i)
		m[i].crop(n);

	for (size_t i = 0; i < this->n; ++i)
	{
		for (size_t j = 0; j < this->n; ++j)
		{
			if (i < n && j < n)
				m[0](i, j) = data[i * this->n + j];
			else if (i < n && j >= n)
				m[1](i, j - n) = data[i * this->n + j];
			else if (i >= n && j < n)
				m[2](i - n, j) = data[i * this->n + j];
			else
				m[3](i - n, j - n) = data[i * this->n + j];
		}
	}
	return m;
}

template <typename T>
SMatrix<T> SMatrix<T>::merge(SMatrix<T> *m)
{
	SMatrix<T> tmp(m[0].n * 2);

	for (size_t i = 0; i < tmp.n; ++i)
	{
		for (size_t j = 0; j < tmp.n; ++j)
		{
			if (i < m[0].n && j < m[0].n)
				tmp(i, j) = m[0](i, j);
			else if (i < m[0].n && j >= m[0].n)
				tmp(i, j) = m[1](i, j - m[0].n);
			else if (i >= m[0].n && j < m[0].n)
				tmp(i, j) = m[2](i - m[0].n, j);
			else
				tmp(i, j) = m[3](i - m[0].n, j - m[0].n);
		}
	}
	return tmp;
}

template <typename T>
void SMatrix<T>::swapRows(size_t i, size_t j)
{
	if (i >= n || j >= n)
		throw std::invalid_argument("Index out of range!");
	if (i == j)
		return;

	for (size_t k = 0; k < n; ++k)
	{
		T tmp = data[i * n + k];
		data[i * n + k] = data[j * n + k];
		data[j * n + k] = tmp;
	}
}

template <typename T>
size_t SMatrix<T>::rank() const
{
	size_t rank = 0;

	SMatrix<T> tmp = gaussEl();

	for (size_t i = 0; i < n; ++i)
		if (tmp(i, i) != 0)
			++rank;
		else break;

	return rank;
}

template <typename T>
SMatrix<T> SMatrix<T>::identity(size_t n)
{
	SMatrix<T> tmp(n);
	tmp.fill(0);
	for (size_t i = 0; i < n; ++i)
		tmp(i, i) = 1;
	return tmp;
}

template <typename T>
T SMatrix<T>::determinant() const
{
	if (n == 1)
		return data[0];
	if (n == 2)
		return data[0] * data[3] - data[1] * data[2];
	// Sarrus rule
	if (n == 3)
		return data[0] * data[4] * data[8] +
			   data[1] * data[5] * data[6] +
			   data[2] * data[3] * data[7] -
			   data[2] * data[4] * data[6] -
			   data[1] * data[3] * data[8] -
			   data[0] * data[5] * data[7];

	T det = 0;
	for (size_t i = 0; i < n; ++i)
	{
		SMatrix<T> tmp(n - 1);
		for (size_t j = 1; j < n; ++j)
		{
			for (size_t k = 0; k < n; ++k)
			{
				if (k < i)
					tmp(j - 1, k) = data[j * n + k];
				else if (k > i)
					tmp(j - 1, k - 1) = data[j * n + k];
			}
		}
		det += data[i] * tmp.determinant() * (i % 2 == 0 ? 1 : -1);
	}
	return det;
}

template <typename T>
SMatrix<T> SMatrix<T>::transpose() const
{
	SMatrix<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			tmp(i, j) = data[j * n + i];
		}
	}
	return tmp;
}

template <typename T>
SMatrix<T> SMatrix<T>::reduce(size_t i, size_t j) const
{
	SMatrix<T> tmp(n - 1);
	for (size_t k = 0; k < n; ++k)
	{
		for (size_t l = 0; l < n; ++l)
		{
			if (k < i && l < j)
				tmp(k, l) = data[k * n + l];
			else if (k < i && l > j)
				tmp(k, l - 1) = data[k * n + l];
			else if (k > i && l < j)
				tmp(k - 1, l) = data[k * n + l];
			else if (k > i && l > j)
				tmp(k - 1, l - 1) = data[k * n + l];
		}
	}
	return tmp;
}

template <typename T>
SMatrix<T> SMatrix<T>::adjoint() const
{
	SMatrix<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			tmp(i, j) = reduce(i, j).determinant() * (i % 2 == j % 2 ? 1 : -1);
		}
	}
	return tmp.transpose();
}

template <typename T>
SMatrix<T> SMatrix<T>::gaussEl() const
{
	SMatrix<T> tmp(n);
	tmp = *this;

	for (size_t i = 0; i < tmp.n; ++i)
	{
		if (tmp(i, i) != 0)
		{
			for (size_t j = i + 1; j < tmp.n; ++j)
			{
				T factor = tmp(j, i) / tmp(i, i);
				for (size_t k = i; k < tmp.n; ++k)
				{
					tmp(j, k) -= factor * tmp(i, k);
				}
			}
		}
		else
		{
			bool non_zero_found = false;
			for (size_t j = i + 1; j < tmp.n; ++j)
			{
				if (tmp(j, i) != 0)
				{
					non_zero_found = true;
					tmp.swapRows(i, j);
					break;
				}
			}
			if (!non_zero_found)
			{
				for (size_t j = i; j < tmp.n; ++j)
				{
					tmp(j, i) = 0;
				}
			}
		}
	}
	return tmp;
}