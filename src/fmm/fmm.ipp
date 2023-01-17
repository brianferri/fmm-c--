template <typename T>
SMatrix<T> fmm::multIt(const SMatrix<T> &A, const SMatrix<T> &B)
{
	if (A.size() != B.size())
		throw std::range_error("Sizes are incompatible! > IT");

	size_t n = A.size();
	SMatrix<T> C(n);
	C.fill(0);

	for (size_t i = 0; i < n; ++i)
	{
		for (size_t k = 0; k < n; ++k)
		{
			for (size_t j = 0; j < n; ++j)
			{
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return C;
}

template <typename T>
SMatrix<T> fmm::multSt(const SMatrix<T> &A, const SMatrix<T> &B)
{
	if (A.size() != B.size())
		throw std::range_error("Sizes are incompatible! > ST");
	if (A.size() == 1)
		return { { A(0, 0) * B(0, 0) } };

	SMatrix<T> *a = A.split();
	SMatrix<T> *b = B.split();

	SMatrix<T> p1 = multSt(a[0] + a[3], b[0] + b[3]);
	SMatrix<T> p2 = multSt(a[2] + a[3], b[0]);
	SMatrix<T> p3 = multSt(a[0], b[1] - b[3]);
	SMatrix<T> p4 = multSt(a[3], b[2] - b[0]);
	SMatrix<T> p5 = multSt(a[0] + a[1], b[3]);
	SMatrix<T> p6 = multSt(a[2] - a[0], b[0] + b[1]);
	SMatrix<T> p7 = multSt(a[1] - a[3], b[2] + b[3]);

	SMatrix<T> *c = new SMatrix<T>[4];
	c[0] = p1 + p4 - p5 + p7;
	c[1] = p3 + p5;
	c[2] = p2 + p4;
	c[3] = p1 - p2 + p3 + p6;

	delete[] a;
	delete[] b;

	return SMatrix<T>::merge(c);
}