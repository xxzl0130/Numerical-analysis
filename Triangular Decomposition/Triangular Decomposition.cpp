#include "Triangular Decomposition.h"
using namespace Eigen;

void DoolittleLU(MatrixXd A, MatrixXd& L, MatrixXd& U)
{
	if(A.rows() != A.cols())
	{
		return;
	}
	const auto n = A.rows();
	L.setIdentity(n, n);
	U.setZero(n, n);
	for(auto k = 0;k < n;++k)
	{
		for(auto j = k;j < n;++j)
		{
			double sum = 0;
			for(auto t = 0;t < k;++t)
			{
				sum += L(k, t) * U(t, j);
			}
			U(k, j) = A(k, j) - sum;
		}
		for(auto i = k + 1;i < n;++i)
		{
			double sum = 0;
			for(auto t = 0;t < k;++t)
			{
				sum += L(i, t) * U(t, k);
			}
			L(i, k) = (A(i, k) - sum) / U(k, k);
		}
	}
}

void CroutLU(MatrixXd A, MatrixXd& L, MatrixXd& U)
{

	if (A.rows() != A.cols())
	{
		return;
	}
	const auto n = A.rows();
	U.setIdentity(n, n);
	L.setZero(n, n);
	for (auto k = 0; k < n; ++k)
	{
		for (auto i = k; i < n; ++i)
		{
			double sum = 0;
			for (auto t = 0; t < k; ++t)
			{
				sum += L(i, t) * U(t, k);
			}
			L(i, k) = A(i, k) - sum;
		}
		for (auto j = k + 1; j < n; ++j)
		{
			double sum = 0;
			for (auto t = 0; t < k; ++t)
			{
				sum += L(k, t) * U(t, j);
			}
			U(k, j) = (A(k, j) - sum) / L(k, k);
		}
	}
}
