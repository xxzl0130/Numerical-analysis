#include "Gaussian Elimination.h"
using namespace Eigen;

VectorXd SequenceGaussian(MatrixXd A)
{
	if (A.cols() - A.rows() != 1)
	{
		return VectorXd();
	}
	const auto n = A.rows();
	VectorXd B = A.col(A.cols() - 1);
	for (auto k = 0; k < n - 1; ++k)
	{
		if (A(k, k) == 0.0)
		{
			return VectorXd();
		}
		for (auto i = k + 1; i < n; ++i)
		{
			auto m = A(i, k) / A(k, k);
			for (auto j = k + 1; j < n; ++j)
			{
				A(i, j) -= m * A(k, j);
			}
			B(i) -= m * B(k);
		}
	}

	VectorXd X;
	X.resize(n);
	X(n - 1) = B(n - 1) / A(n - 1, n - 1);
	for(auto i = n - 2;i >= 0;--i)
	{
		double sum = 0;
		for(auto j = i + 1;j < n;++j)
		{
			sum += A(i, j) * X(j);
		}
		X(i) = (B(i) - sum) / A(i, i);
	}

	return X;
}

VectorXd mainElementGaussian(MatrixXd A)
{
	if (A.cols() - A.rows() != 1)
	{
		return VectorXd();
	}
	const auto n = A.rows();
	VectorXd B = A.col(A.cols() - 1);

	for (auto k = 0; k < n; ++k)
	{
		int index = k;
		double max = abs(A(index, k));
		for (auto i = k + 1; i < n; ++i)
		{
			if (abs(A(i, k)) > max)
			{
				max = abs(A(i, k));
				index = i;
			}
		}
		A.row(k).swap(A.row(index));
		B.row(k).swap(B.row(index));
		for (auto i = k + 1; i < n; ++i)
		{
			auto m = A(i, k) / A(k, k);
			for (auto j = k + 1; j < n; ++j)
			{
				A(i, j) -= m * A(k, j);
			}
			B(i) -= m * B(k);
		}
	}

	VectorXd X;
	X.resize(n);
	X(n - 1) = B(n - 1) / A(n - 1, n - 1);
	for (auto i = n - 2; i >= 0; --i)
	{
		double sum = 0;
		for (auto j = i + 1; j < n; ++j)
		{
			sum += A(i, j) * X(j);
		}
		X(i) = (B(i) - sum) / A(i, i);
	}

	return X;
}
