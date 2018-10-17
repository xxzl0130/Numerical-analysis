#include "Triangular Decomposition.h"
#include <algorithm>
#include <iostream>
using namespace Eigen;
using std::swap;

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

void MainElementDoolittleLU(MatrixXd A, MatrixXd& L, MatrixXd& U)
{
	if (A.rows() != A.cols())
	{
		return;
	}
	const auto n = A.rows();
	L.setIdentity(n, n);
	U.setZero(n, n);
	VectorXi M(n);
	VectorXd S(n);
	for(auto i = 0;i < n;++i)
	{
		M(i) = i;
	}
	for(auto k = 0;k < n;++k)
	{
		for(auto i = k;i < n;++i)
		{
			double sum = 0;
			for(auto t = 0;t < k;++t)
			{
				sum += L(i, t) * U(t, k);
			}
			S(i) = A(i, k) - sum;
		}
		auto max = abs(S(k));
		auto index = k;
		for(auto i = k + 1;i < n;++i)
		{
			if(abs(S(i)) > max)
			{
				max = abs(S(i));
				index = i;
			}
		}
		//M(k) = index;
		if(index != k)
		{
			L.row(k).swap(L.row(index));
			A.row(k).swap(A.row(index));
			/*for(auto t = 0;t < k;++t)
			{
				swap(L(k, t), L(index, t));
			}
			for(auto t = k;t < n;++t)
			{
				swap(A(k, t), A(index, t));
			}*/
			//S.row(k).swap(S.row(index));
			swap(S(k), S(index));
			M.row(k).swap(M.row(index));
		}
		
		U(k, k) = S(k);
		for (auto j = k; j < n; ++j)
		{
			double sum = 0;
			for (auto t = 0; t < k; ++t)
			{
				sum += L(k, t) * U(t, j);
			}
			U(k, j) = A(k, j) - sum;
		}
		for (auto i = k + 1; i < n; ++i)
		{
			L(i, k) = S(i) / U(k, k);
		}
	}
	std::cout << M << std::endl;
	MatrixXd _L, _U;
	_L.setZero(n, n);
	_U.setZero(n, n);
	for(auto i = 0;i < n;++i)
	{
		_L.row(M(i)) = L.row(i);
		_U.row(M(i)) = U.row(i);
	}
	//L = _L;
	//U = _U;
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
