template <typename T>
SquareMatrix<T> fmm::mult(const SquareMatrix<T> &A, const SquareMatrix<T> &B)
{
	if (A.size() != B.size())
		throw std::range_error("Sizes are incompatible! > IT");
	size_t n = A.size();
	SquareMatrix<T> C(n);
	C.fill(0);
	for (size_t i = 0; i < n; ++i)
		for (size_t k = 0; k < n; ++k)
			for (size_t j = 0; j < n; ++j)
				C(i, j) += A(i, k) * B(k, j);
	return C;
}

template <typename T>
SquareMatrix<T> fmm::strassen(const SquareMatrix<T> &A, const SquareMatrix<T> &B)
{
	if (A.size() != B.size())
		throw std::range_error("Sizes are incompatible! > ST");
	if (A.size() == 2)
	{
		T m1 = (A(0, 0) + A(1, 1)) * (B(0, 0) + B(1, 1));
		T m2 = (A(1, 0) + A(1, 1)) * B(0, 0);
		T m3 = A(0, 0) * (B(0, 1) - B(1, 1));
		T m4 = A(1, 1) * (B(1, 0) - B(0, 0));
		T m5 = (A(0, 0) + A(0, 1)) * B(1, 1);
		T m6 = (A(1, 0) - A(0, 0)) * (B(0, 0) + B(0, 1));
		T m7 = (A(0, 1) - A(1, 1)) * (B(1, 0) + B(1, 1));
		return { { m1 + m4 - m5 + m7, m3 + m5 },
				 { m2 + m4, m1 - m2 + m3 + m6 } };
	} /* Just in case */ else if (A.size() == 1) {
		return { { A(0, 0) * B(0, 0) } };
	} // Base case

	SquareMatrix<T> *a = A.split();
	SquareMatrix<T> *b = B.split();

	SquareMatrix<T> p1 = strassen(a[0] + a[3], b[0] + b[3]);
	SquareMatrix<T> p2 = strassen(a[2] + a[3], b[0]);
	SquareMatrix<T> p3 = strassen(a[0], b[1] - b[3]);
	SquareMatrix<T> p4 = strassen(a[3], b[2] - b[0]);
	SquareMatrix<T> p5 = strassen(a[0] + a[1], b[3]);
	SquareMatrix<T> p6 = strassen(a[2] - a[0], b[0] + b[1]);
	SquareMatrix<T> p7 = strassen(a[1] - a[3], b[2] + b[3]);

	SquareMatrix<T> *c = new SquareMatrix<T>[4];
	c[0] = p1 + p4 - p5 + p7;
	c[1] = p3 + p5;
	c[2] = p2 + p4;
	c[3] = p1 - p2 + p3 + p6;

	delete[] a;
	delete[] b;

	SquareMatrix<T> m = SquareMatrix<T>::merge(c);
	delete[] c;

	m.trim(A.size());
	return m;
}