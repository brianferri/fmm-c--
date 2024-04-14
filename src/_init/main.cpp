#include "../../include/fmm/fmm.tpp"

using namespace fmm;

int main()
{
	/* Init */
	SMatrix<double> A({{1, 3, 2},
					   {2, 0, 1},
					   {1, 2, 0}});

	SMatrix<double> B({{3, 0, 1},
					   {1, 0, 1},
					   {0, 1, 2}});

	/* Test Section */
	try
	{
		std::cout << "-----------Start-----------" << std::endl;
		std::cout << A << std::endl;
		std::cout << B << std::endl;

		std::cout << "-----------Crop-----------" << std::endl;
		A.crop(4);
		std::cout << A << std::endl;

		std::cout << "-----------Trim-----------" << std::endl;
		A.trim(3);
		std::cout << A << std::endl;

		std::cout << "-----------Split-----------" << std::endl;
		SMatrix<double> *C = B.split();
		for (size_t i = 0; i < 4; ++i)
			std::cout << C[i] << std::endl;

		std::cout << "-----------Assign-----------" << std::endl;
		for (size_t i = 0; i < 4; ++i)
		{
			C[i] = A;
			std::cout << C[i] << std::endl;
		}

		std::cout << "-----------Merge-----------" << std::endl;
		SMatrix<double> D(4);
		D = SMatrix<double>::merge(C);
		std::cout << D << std::endl;

		std::cout << "-----------SwapRows-----------" << std::endl;
		D.swapRows(0, 1);
		std::cout << D << std::endl;

		std::cout << "-----------Sum-----------" << std::endl;
		SMatrix<double> E = C[0] + C[1];
		std::cout << E << std::endl;

		std::cout << "-----------Sub-----------" << std::endl;
		E = C[0] - C[1];
		std::cout << E << std::endl;

		std::cout << "-----------MultbyScalar-----------" << std::endl;
		E = C[0] * 3;
		E = 3 * C[0];
		std::cout << E << std::endl;

		std::cout << "-----------Determinant-----------" << std::endl;
		std::cout << E.determinant() << std::endl;

		std::cout << "-----------Transpose-----------" << std::endl;
		std::cout << B.transpose() << std::endl;

		std::cout << "-----------Minor-----------" << std::endl;
		std::cout << B.minor(1, 1) << std::endl;

		std::cout << "-----------Adjoint-----------" << std::endl;
		std::cout << B.adjoint() << std::endl;

		std::cout << "-----------Row Echelon Form-----------" << std::endl;
		std::cout << A.REF() << std::endl;

		std::cout << "-----------Rank-----------" << std::endl;
		std::cout << A.rank() << std::endl;
		std::cout << B.rank() << std::endl;
		std::cout << (A + B).rank() << std::endl;

		std::cout << "-----------Naive-----------" << std::endl;
		D = multIt(A, B);
		std::cout << D << std::endl;

		std::cout << "-----------Strassen-----------" << std::endl;
		E = multSt(A, B);
		E.trim(A.size()); // TODO: Fix this (trim() should be called automatically)
		std::cout << E << std::endl;

		std::cout << "-----------3 Dimensional Matrix-----------" << std::endl;
		SMatrix<SMatrix<double> > F(2);
		F[0] = A;
		F[1] = B;
		std::cout << F << std::endl;

		std::cout << "-----------3 Dimensional Matrix By List Initializer-----------" << std::endl;
		SMatrix<SMatrix<double> > G({{A, B},
									 {B, A}});
		std::cout << G << std::endl;

		std::cout << "-----------3 Dimensional Matrix By Default Constructor-----------" << std::endl;
		SMatrix<SMatrix<double> > H = SMatrix<SMatrix<double> >(2);
		H[0] = A;
		H[1] = B;
		std::cout << H << std::endl;

		std::cout << "-----------3 Dimensional Matrix By Pointer-----------" << std::endl;
		SMatrix<SMatrix<double> >* I = new SMatrix<SMatrix<double> >(2);
		(*I)[0] = A;
		(*I)[1] = B;
		std::cout << *I << std::endl;

		std::cout << "-----------End-----------" << std::endl;
	}
	catch (std::runtime_error err)
	{
		std::cout << "Error: " << err.what() << std::endl;
	}
}