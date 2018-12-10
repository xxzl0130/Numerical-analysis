#include <iostream>
#include <cmath>
#include <cstdio>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace Eigen;

#define EPS			1e-12	//精度
#define MAX_LOOP	1000	//最大迭代次数
#define DIM			4		//方程维数
#define N_GIVEN_T	6		//t的个数
#define H_T			0.2		//t的间隔
#define N_GIVEN_U	6		//u的个数
#define H_U			0.4		//u的间隔
#define N_X			11		//xi的个数
#define H_X			0.08	//xi的间隔
#define N_Y			21		//yi的个数
#define H_Y			0.05	//yi的间隔
typedef Matrix<double, 6, 1> Vector6d;

double givenZ[N_GIVEN_T][N_GIVEN_U] = {
	{-0.5, -0.34, 0.14, 0.94, 2.06, 3.5},
	{-0.42,-0.5, -0.26, 0.3,  1.18, 2.38},
	{-0.18,-0.5, -0.5, -0.18, 0.46, 1.42},
	{ 0.22,-0.34,-0.58,-0.5, -0.1,  0.62},
	{ 0.78,-0.02,-0.5, -0.66,-0.5, -0.02},
	{ 1.5,  0.46,-0.26,-0.66,-0.74,-0.5}
};
double givenT[N_GIVEN_T] = { 0,0.2,0.4,0.6,0.8,1.0 };
double givenU[N_GIVEN_U] = { 0,0.4,0.8,1.2,1.6,2.0 };
double x[N_X], y[N_Y], fxy[N_X][N_Y];

/**
 * \brief 非线性方程函数值
 * \param X 变量组[x,y,t,u,v,w]
 * \return 4个函数值
 */
Vector4d F(Vector6d X);
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
double lagrangeInterpolate(double t, double u);
/**
 * \brief 二元曲面拟合
 * \param x 自变量x
 * \param y 自变量y
 * \param fxy 函数值
 * \param k 拟合次数
 * \return 系数矩阵
 */
MatrixXd surfacesFit(const double x[N_X], const double y[N_Y], const double fxy[N_X][N_Y],int k);
double power(double x, int a);
double testError(const MatrixXd& c, const double x[N_X], const double y[N_Y], const double fxy[N_X][N_Y]);
/**
 * \brief p(x,y)拟合函数
 * \param c 系数矩阵
 * \param x x
 * \param y y
 * \return p(x,y)
 */
double p(const MatrixXd& c, double x, double y);

int main()
{
	fstream fout;
	fout.open("output.txt", ios_base::out);

	for(auto i = 0;i < N_X;++i)
		x[i] = H_X * i;
	for (auto j = 0; j < N_Y; ++j)
		y[j] = 0.5 + H_Y * j;

	Vector4d initState = Vector4d::Ones();
	fout << "x,y,f" << endl;
	for(auto i = 0;i < N_X;++i)
	{
		for(auto j = 0;j < N_Y;++j)
		{
			auto var = newton(initState, x[i], y[j]);
			fxy[i][j] = lagrangeInterpolate(var[0], var[1]);
			fout << setprecision(0) << resetiosflags(ios::scientific);
			fout << x[i] << "," << y[j] << ",";
			fout << setprecision(12) << setiosflags(ios::scientific);
			fout << fxy[i][j] << endl;
		}
	}

	fout << setprecision(12) << setiosflags(ios::scientific);
	fout << endl << "k,delta\n";
	MatrixXd Crs;
	for(auto k = 0;;++k)
	{
		Crs = surfacesFit(x, y, fxy, k);
		auto err = testError(Crs, x, y, fxy);
		fout << setprecision(0) << resetiosflags(ios::scientific);
		fout << k << ",";
		fout << setprecision(12) << setiosflags(ios::scientific);
		fout << err << endl;
		if (err < 1e-7)
			break;
	}
	fout << endl << "Crs:\n" << Crs << endl;

	fout << endl << "x*,y*,f,p" << endl;
	for(auto i = 1;i <= 8;++i)
	{
		auto x = 0.1 * i;
		for(auto j = 1;j <= 5;++j)
		{
			auto y = 0.5 + 0.2 * j;
			auto var = newton(initState, x, y);
			auto fxy = lagrangeInterpolate(var[0], var[1]);
			auto pxy = p(Crs, x, y);
			fout << setprecision(0) << resetiosflags(ios::scientific);
			fout << x << "," << y << ",";
			fout << setprecision(12) << setiosflags(ios::scientific);
			fout << fxy << "," << pxy << endl;
		}
	}
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

double lagrangeInterpolate(double t, double u)
{
	double p = 0.0;
	int i = 1, j = 1;

	if(t <= givenT[1] + H_T / 2)
	{
		i = 1;
	}
	else if(t > givenT[N_GIVEN_T - 2] - H_T / 2)
	{
		i = N_GIVEN_T - 2;
	}
	else
	{
		for(auto k = 2;k < N_GIVEN_T - 2;++k)
		{
			if(givenT[k] - H_T / 2 < t && t <= givenT[k] + H_T / 2)
			{
				i = k;
			}
		}
	}

	if (u <= givenU[1] + H_U / 2)
	{
		j = 1;
	}
	else if (u > givenU[N_GIVEN_U - 2] - H_U / 2)
	{
		j = N_GIVEN_U - 2;
	}
	else
	{
		for (auto k = 2; k < N_GIVEN_U - 2; ++k)
		{
			if (givenU[k] - H_U / 2 < u && u <= givenU[k] + H_U / 2)
			{
				j = k;
			}
		}
	}

	for(auto k = i - 1;k <= i + 1;++k)
	{
		for(auto r = j - 1;r <= j + 1;++r)
		{
			auto tmp = 1.0;
			for(auto s = i -1;s <= i + 1;++s)
			{
				if(s == k)
					continue;
				tmp *= (t - givenT[s]) / (givenT[k] - givenT[s]);
			}
			for(auto s = j - 1;s <= j + 1;++s)
			{
				if(s == r)
					continue;
				tmp *= (u - givenU[s]) / (givenU[r] - givenU[s]);
			}
			tmp *= givenZ[k][r];
			p += tmp;
		}
	}

	return p;
}

MatrixXd surfacesFit(const double x[11], const double y[21], const double fxy[11][21], int k)
{
	MatrixXd Crs;
	MatrixXd f = Eigen::Map<Matrix<double, N_X, N_Y, RowMajor>>((double*)fxy);
	//构造Bij = phi_j(xi)
	MatrixXd B;
	B.resize(N_X, k + 1);
	for (auto i = 0; i < N_X; ++i)
		for (auto j = 0; j <= k; ++j)
			B(i, j) = pow(x[i], j);
	MatrixXd BTB = B.transpose() * B;
	MatrixXd alpha;
	alpha.resize(k + 1, N_Y);
	for (auto j = 0; j < N_Y; ++j)
	{
		VectorXd u = f.col(j);
		alpha.col(j) = MainElementGaussian(BTB, B.transpose() * u);
	}

	//构造Gij = phi_j(yi)
	MatrixXd G;
	G.resize(N_Y, k + 1);
	for (auto i = 0; i < N_Y; ++i)
		for (auto j = 0; j <= k; ++j)
			G(i, j) = pow(y[i], j);
	MatrixXd GTG = G.transpose() * G;

	MatrixXd beta;
	beta.resize(k + 1, N_Y);
	for (auto j = 0; j < N_Y; ++j)
	{
		VectorXd q = alpha.row(0);
		beta.col(j) = MainElementGaussian(GTG, G.row(j));
	}

	Crs = alpha * beta.transpose();
	return Crs;
}

double power(double x, int a)
{
	double s = 1.0;
	while(a)
	{
		if (a & 1)
			s *= x;
		x *= x;
		a >>= 1;
	}
	return s;
}

double testError(const MatrixXd& c, const double x[11], const double y[21], const double fxy[11][21])
{
	double err = 0.0;
	auto k = c.rows();
	for(auto i = 0;i < N_X;++i)
	{
		for(auto j = 0;j < N_Y;++j)
		{
			err += power(fxy[i][j] - p(c, x[i], y[j]), 2);
		}
	}
	return err;
}

double p(const MatrixXd& c, double x, double y)
{
	auto p = 0.0;
	auto k = c.rows();
	for (auto r = 0; r < k; ++r)
	{
		for (auto s = 0; s < k; ++s)
		{
			p += c(r, s) * power(x, r) * power(y, s);
		}
	}
	return p;
}
