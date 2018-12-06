#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

#define EPS			1e-12	//精度
#define DIM			10		//维度
#define MAX_LOOP	5000	//最大迭代次数

MatrixXd A, Q, R;

//初始化A矩阵
void init(MatrixXd& M);
/**
 * \brief 拟上三角化
 * \param A 输入矩阵
 * \return 迭代到n-1次的矩阵
 */
MatrixXd quasiTriangular(MatrixXd A);
//浮点数判0
inline bool realZero(double x)
{
	return abs(x) <= EPS;
}
//符号函数
inline double sgn(double x)
{
	return x > 0 ? 1 : (x < 0 ? -1 : 0);
}
void checkZero(MatrixXd& M);
MatrixXd checkZero(MatrixXd M);
/**
 * \brief 带双步位移的QR分解
 * \param A 输入矩阵
 * \return 特征值
 */
MatrixXd twoStepQrTrans(MatrixXd A);
/**
 * \brief 二阶矩阵特征根
 * \param A 输入矩阵
 * \return 特征根，每行一个
 */
Matrix2d twoOrderLambda(Matrix2d A);
/**
 * \brief QR分解
 * \param A 输入矩阵
 * \param M 输入矩阵
 * \return A(k+1) = Qk^T*Ak*Qk
 */
MatrixXd QR(MatrixXd M, MatrixXd A);
/**
 * \brief QR分解
 * \param A 输入矩阵
 * \param Q 输出Q
 * \param R 输出R
 */
void QR(MatrixXd A, MatrixXd &Q, MatrixXd& R);
/**
 * \brief 求特征向量
 * \param A 输入矩阵
 * \param lambda 特征值
 * \return 特征向量
 */
VectorXd eigenVector(MatrixXd A, double lambda);
/**
 * \brief 选主元的高斯消元
 * \param A 增广矩阵
 * \return 解
 */
VectorXd mainElementGaussian(MatrixXd A);

int main()
{
	cout << scientific << setprecision(12);//12位科学计数法
	
	init(A);
	A = checkZero(quasiTriangular(A));
	cout << "拟上三角阵：\n" << A << endl;
	
	QR(A, Q, R);
	cout << "Q矩阵为：\n" << Q << endl;
	cout << "R矩阵为：\n" << R << endl;

	auto lambda = twoStepQrTrans(A);
	cout << "特征值为：\n" << lambda << endl;

	for(auto i = 0;i < lambda.rows();++i)
	{
		if(realZero(lambda(i,1)))
		{//实根
			cout << lambda(i, 0) << "对应的特征向量为：\n" 
			<< eigenVector(A, lambda(i, 0)).transpose() << endl;
		}
	}

	system("pause");
}

void init(MatrixXd& M)
{
	M.resize(DIM, DIM);
	for(auto i = 1;i <= DIM;++i)
	{
		for(auto j = 1;j <= DIM;++j)
		{
			if(i == j)
			{
				M(i - 1,j - 1) = 1.52 * cos(i + 1.2 * j);
			}
			else
			{
				M(i - 1, j - 1) = sin(0.5 * i + 0.2 * j);
			}
		}
	}
}

MatrixXd quasiTriangular(MatrixXd A)
{
	double d, c, h;
	VectorXd u(DIM);

	for(auto r = 0;r < DIM - 2;++r)
	{
		auto allZero = true;
		for(auto i = r + 2;i < DIM;++i)
		{
			allZero = allZero && realZero(A(i,r));
		}
		if(!allZero)
		{
			d = 0.0;
			for(auto i = r + 1;i < DIM;++i)
			{
				d += A(i, r) * A(i, r);
			}
			d = sqrt(d);
			if(realZero(A(r + 1, r)))
			{
				c = d;
			}
			else
			{
				c = -sgn(A(r + 1,r)) * d;
			}
			h = c * c - c * A(r + 1, r);
			u.setZero();
			u[r + 1] = A(r + 1, r) - c;
			for(auto i = r + 2;i < DIM;++i)
			{
				u[i] = A(i, r);
			}
			VectorXd p = A.transpose() * u / h;
			VectorXd q = A * u / h;
			double t = p.dot(u) / h;
			VectorXd w = q - t * u;
			A = A - w * u.transpose() - u * p.transpose();
		}
	}
	return A;
}

void checkZero(MatrixXd& M)
{
	for(auto i = 0;i < M.rows();++i)
	{
		for(auto j = 0;j < M.cols();++j)
		{
			if(realZero(M(i,j)))
			{
				M(i, j) = 0.0;
			}
		}
	}
}

MatrixXd checkZero(MatrixXd M)
{
	for (auto i = 0; i < M.rows(); ++i)
	{
		for (auto j = 0; j < M.cols(); ++j)
		{
			if (realZero(M(i, j)))
			{
				M(i, j) = 0.0;
			}
		}
	}
	return M;
}

MatrixXd twoStepQrTrans(MatrixXd A)
{
	auto n = A.rows();
	auto m = n - 1;
	auto k = 1,i = 0;
	MatrixXd lambda = MatrixXd::Zero(A.rows(), 2);

	while(k < MAX_LOOP)
	{
		//(3)
		if(abs(A(m,m - 1)) <= EPS)//得到一个特征值
		{
			lambda(i, 0) = A(m, m);
			lambda(i, 1) = 0;
			--m;
			A = A.block(0, 0, m + 1, m + 1);//Ak=[aij]mxm
			++i;
		}
		else
		{
			goto step5;
		}
		//(4)
		step4:
		if (m == 0)
		{
			lambda(i, 0) = A(m, m);
			lambda(i, 1) = 0;
			break;//转11
		}
		else if (m < 0)
		{
			break;//转11
		}
		else
		{
			continue;//转3
		}
		step5:
		//(5)
		auto l = twoOrderLambda(A.block(m - 1, m - 1, 2, 2));
		//(6)
		if(m == 1)
		{
			lambda.row(i++) = l.row(0);
			lambda.row(i++) = l.row(1);
			break;//转11
		}
		//(7)
		if(abs(A(m - 1,m - 2)) <= EPS)
		{
			lambda.row(i++) = l.row(0);
			lambda.row(i++) = l.row(1);
			m -= 2;
			A = A.block(0, 0, m + 1, m + 1);//Ak=[aij]mxm
			++k;
			goto step4;
		}
		else if(k >= MAX_LOOP)//(8)
		{
			break;
		}
		else
		{
			//(9)
			auto s = A(m - 1, m - 1) + A(m, m);
			auto t = A(m - 1, m - 1) * A(m, m) - A(m, m - 1)*A(m - 1, m);
			MatrixXd M = A * A - s * A + t * MatrixXd::Identity(m + 1, m + 1);
			A = QR(M, A);
			++k;
		}
	}
	cout << "QR迭代后的A阵为：" << A << endl;
	return lambda;
}

Matrix2d twoOrderLambda(Matrix2d A)
{
	Matrix2d ans = Matrix2d::Zero();

	auto b = -(A(0, 0) + A(1, 1));
	auto c = A.determinant();//行列式
	auto det = b * b - 4 * c;
	if(det >= 0)
	{//2实根
		ans(0, 0) = (-b + sqrt(det)) / 2;
		ans(1, 0) = (-b - sqrt(det)) / 2;
	}
	else
	{//2共轭虚根
		ans(0, 0) = ans(1, 0) = -b / 2;
		ans(0, 1) = sqrt(abs(det)) / 2;
		ans(1, 1) = -sqrt(abs(det)) / 2;
	}

	return ans;
}

MatrixXd QR(MatrixXd M, MatrixXd A)
{
	auto &B = M;
	auto &C = A;
	auto m = A.rows();
	for(auto r = 0;r < m - 1;++r)
	{
		bool allZero = true;
		for(auto i = r + 1;i < m;++i)
		{
			allZero = allZero && realZero(B(i, r));
		}
		if(allZero)
		{
			continue;
		}
		auto tmp = 0.0;
		for(auto i = r;i < m;++i)
		{
			tmp += B(i, r) * B(i, r);
		}
		auto d = sqrt(tmp);
		auto c = -sgn(B(r, r)) * d;
		if(realZero(B(r,r)))
		{
			c = d;
		}
		auto h = c * c - c * B(r, r);
		VectorXd u = VectorXd::Zero(m);
		u(r) = B(r, r) - c;
		for(auto i = r + 1;i < m;++i)
		{
			u(i) = B(i, r);
		}
		VectorXd v = B.transpose() * u / h;
		B = B - u * v.transpose();
		VectorXd p = C.transpose() * u / h;
		VectorXd q = C * u / h;
		double t = p.dot(u) / h;
		VectorXd w = q - t * u;
		C = C - w * u.transpose() - u * p.transpose();
	}
	return C;
}

void QR(MatrixXd A, MatrixXd& Q, MatrixXd& R)
{
	auto m = A.rows();
	Q = MatrixXd::Identity(m, m);
	for (auto r = 0; r < m - 1; ++r)
	{
		bool allZero = true;
		for (auto i = r + 1; i < m; ++i)
		{
			allZero = allZero && realZero(A(i, r));
		}
		if (allZero)
		{
			continue;
		}
		auto tmp = 0.0;
		for (auto i = r; i < m; ++i)
		{
			tmp += A(i, r) * A(i, r);
		}
		auto d = sqrt(tmp);
		auto c = -sgn(A(r, r)) * d;
		if (realZero(A(r, r)))
		{
			c = d;
		}
		auto h = c * c - c * A(r, r);
		VectorXd u = VectorXd::Zero(m);
		u(r) = A(r, r) - c;
		for (auto i = r + 1; i < m; ++i)
		{
			u(i) = A(i, r);
		}
		VectorXd w = Q * u;
		Q = Q - w * u.transpose() / h;
		VectorXd p = A.transpose() * u / h;
		A = A - u * p.transpose();
	}
	R = A;
}

VectorXd eigenVector(MatrixXd A, double lambda)
{
	A = A - lambda * MatrixXd::Identity(A.rows(), A.cols());
	auto B = VectorXd::Zero(A.rows());
	MatrixXd M;
	M.resize(A.rows(), A.rows() + 1);
	M.block(0, 0, A.rows(), A.rows()) = A;
	M.col(M.cols() - 1) = B;
	return mainElementGaussian(M);
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
	if (abs(A(n - 1, n - 1)) < 1e-8)
	{
		X(n - 1) = 1;
	}
	else
	{
		X(n - 1) = B(n - 1) / A(n - 1, n - 1);
	}
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
