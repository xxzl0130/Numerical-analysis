#include <iostream>
#include <cmath>
#include <cstdio>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#define EPS			1e-12	//精度
#define MAX_LOOP	1000	//最大迭代次数
#define DIM			4		//方程维数
typedef Matrix<double, 6, 1> Vector6d;

/**
 * \brief 非线性方程函数值
 * \param X 变量组[x,y,t,u,v,w]
 * \return 4个函数值
 */
Vector4d F(Vector6d X);
/**
 * \brief 非线性方程函数的导数值
 * \param x 变量组[x,y,t,u,v,w]
 * \return 4个函数值
 */
Vector4d Fd(Vector6d X);
/**
 * \brief 非线性方程参数阵
 * \param X 变量组[t,u,v,w]
 * \return 参数阵
 */
Matrix4d A(Vector4d X);
/**
 * \brief 非线性方程的导数阵
 * \param X 变量组[t,u,v,w]
 * \return 参数阵
 */
Matrix4d Ad(Vector4d X);
/**
 * \brief 牛顿迭代法求解非线性方程组
 * \param var 4个待求参数[t,u,v,w]
 * \param x 已知x
 * \param y 已知y
 * \return 4个解出参数[t,u,v,w]
 */
Vector4d newton(Vector4d var, double x, double y);
/**
 * \brief 选主元的高斯消元法
 * \param A 系数矩阵
 * \param B 值向量
 * \return 解
 */
Eigen::VectorXd MainElementGaussian(Eigen::MatrixXd A,VectorXd B);


int main()
{
	Vector4d var = Vector4d::Ones();
	Vector6d vars;
	cout << (var = newton(var, 0.08, 0.5).transpose()) << endl;
	vars << 0.08, 0.5, var;
	cout << F(vars).transpose() << endl;

	system("pause");
}

Vector4d F(Vector6d X)
{
	Vector4d r;
	auto &x = X[0], &y = X[1], &t = X[2], &u = X[3], &v = X[4], &w = X[5];
	r <<
		0.5 * cos(t) + u + v + w - x - 2.67,
		t + 0.5 * sin(u) + v + w - y - 1.07,
		0.5 * t + u + cos(v) + w - x - 3.74,
		t + 0.5 * u + v + sin(w) - y - 0.79;
	return r;
}

Vector4d Fd(Vector6d X)
{
	Vector4d r;
	auto &x = X[0], &y = X[1], &t = X[2], &u = X[3], &v = X[4], &w = X[5];
	r <<
		- (-0.5 * sin(t) + 1 + 1 + 1 - 1),
		- (1 + 0.5 * cos(u) + 1 + 1 - 1),
		- (0.5 + 1 - sin(v) + 1 - 1),
		- (1 + 0.5 + 1 + cos(w) - 1);
	return r;
}

Matrix4d A(Vector4d X)
{
	Matrix4d m;
	auto &t = X[0], &u = X[1], &v = X[2], &w = X[3];
	m << 0.5 * cos(t), u, v, w,
		 t, 0.5 * sin(u), v, w,
		 0.5 * t, u, cos(v), w,
		 t, 0.5 * u, v, sin(w);
	return m;
}

Matrix4d Ad(Vector4d X)
{
	Matrix4d m = Matrix4d::Ones();
	auto &t = X[0], &u = X[1], &v = X[2], &w = X[3];
	m(0, 0) = -0.5 * sin(t);
	m(1, 1) = 0.5 * cos(u);
	m(2, 0) = 0.5;
	m(2, 2) = -sin(v);
	m(3, 1) = 0.5;
	m(3, 3) = cos(w);
	return m;
}

Vector4d newton(Vector4d var, double x, double y)
{
	for(auto k = 0;k < MAX_LOOP;++k)
	{
		Vector6d vars;
		vars << x, y, var;
		auto temp = MainElementGaussian(Ad(var), -F(vars));//式4.26

		for(auto i = 0;i < 4;++i)
		{
			if(abs(temp[i] / var[i]) < EPS)
			{
				cout << k << endl;
				return var;
			}
		}
		var += temp;
	}

	return var;
}

VectorXd MainElementGaussian(MatrixXd A, VectorXd B)
{
	const auto n = A.rows();

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
