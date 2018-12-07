#include <iostream>   //声明程序所需要输入输出操作的有关信息
#include <iomanip>    //操作符头文件
#include <fstream>
#include <math.h>
using namespace std;  //针对命名空间的指令

const int I = 11;//定义插值数i
const int J = 21;//定义插值数j
const int L = 4;//非线性方程组维数
const int Max = 10;//最大迭代次数
const double precision = 1e-12;

ofstream outfile;//创建输出流对象

void initial(double f[][6]) //初始化数表
{
	f[0][0] = -0.5; f[0][1] = -0.34; f[0][2] = 0.14; f[0][3] = 0.94; f[0][4] = 2.06; f[0][5] = 3.5;
	f[1][0] = -0.42; f[1][1] = -0.5; f[1][2] = -0.26; f[1][3] = 0.3; f[1][4] = 1.18; f[1][5] = 2.38;
	f[2][0] = -0.18; f[2][1] = -0.5; f[2][2] = -0.5; f[2][3] = -0.18; f[2][4] = 0.46; f[2][5] = 1.42;
	f[3][0] = 0.22; f[3][1] = -0.34; f[3][2] = -0.58; f[3][3] = -0.5; f[3][4] = -0.1; f[3][5] = 0.62;
	f[4][0] = 0.78; f[4][1] = -0.02; f[4][2] = -0.5; f[4][3] = -0.66; f[4][4] = -0.5; f[4][5] = -0.02;
	f[5][0] = 1.5; f[5][1] = 0.46; f[5][2] = -0.26; f[5][3] = -0.66; f[5][4] = -0.74; f[5][5] = -0.5;
}

void Guass_choose(double C[][I], double d[], double x[], int dim) //列主元Guass消去法求方程组的解
{
	double max = 0.0, temp = 0.0;
	double A[I][I] = { 0.0 }, b[I] = { 0.0 };
	int i, j, k;

	for (i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			A[i][j] = C[i][j];
	for (i = 0; i < dim; i++) b[i] = d[i];

	for (k = 0; k < dim - 1; k++)
	{
		max = fabs(A[k][k]);
		for (i = k + 1, j = k; i < dim; i++)
		{
			if (max < fabs(A[i][k]))
			{
				max = fabs(A[i][k]); j += 1;
			}
		}

		if (j != k)
		{
			for (i = 0; i < dim; i++)
			{
				temp = A[j][i]; A[j][i] = A[k][i]; A[k][i] = temp;
			}
			temp = b[j]; b[j] = b[k]; b[k] = temp;
		}

		for (i = k + 1; i < dim; i++)
		{
			for (j = k + 1; j < dim; j++)
				A[i][j] -= A[k][j] * (A[i][k] / A[k][k]);
			b[i] -= b[k] * (A[i][k] / A[k][k]);
		}

	}

	x[dim - 1] = b[dim - 1] / A[dim - 1][dim - 1];
	for (k = dim - 2; k >= 0; k--)
	{
		temp = 0.0;
		for (j = k + 1; j < dim; j++) temp += A[k][j] * x[j];
		x[k] = (b[k] - temp) / A[k][k];
	}
}

void Newton(double variable[], double x, double y) //牛顿法求解非线性方程组
{
	double A[I][I] = { 0.0 }, b[L] = { 0.0 }, temp[L] = { 0.0 };
	double t, u, v, w;
	int k, i, j;
	t = variable[0]; u = variable[1]; v = variable[2]; w = variable[3];
	for (k = 0; k <= Max; k++)
	{
		b[0] = -0.5*cos(t) - u - v - w + x + 2.67;
		b[1] = -t - 0.5*sin(u) - v - w + y + 1.07;
		b[2] = -0.5*t - u - cos(v) - w + x + 3.74;
		b[3] = -t - 0.5*u - v - sin(w) + y + 0.79;
		for (i = 0; i < L; i++)
			for (j = 0; j < L; j++) A[i][j] = 1.0;
		A[0][0] = -0.5*sin(t); A[1][1] = 0.5*cos(u);
		A[2][0] = 0.5; A[2][2] = -sin(v); A[3][1] = 0.5; A[3][3] = cos(w);

		Guass_choose(A, b, temp, L);
		if (fabs(temp[0] / t) < precision || fabs(temp[1] / u) < precision
			|| fabs(temp[2] / v) < precision || fabs(temp[3] / w) < precision)
			break;
		t += temp[0]; u += temp[1]; v += temp[2]; w += temp[3];
	}
	variable[0] = t; variable[1] = u; variable[2] = v; variable[3] = w;
}

double Lagrange(double tx, double uy, double f[][6])//二元分片二次插值
{
	double x[6], y[6];//n=m=5
	int i, j, k, r, t, kx, ry;
	double temp1, temp2, temp3, p22;
	for (i = 0; i < 6; i++)
	{
		y[i] = 0.4*i; x[i] = 0.2*i;//h=0.4,tao=0.2;
	}

	for (i = 2; i < 4; i++)
	{
		if ((tx > (x[i] - 0.1)) && (tx <= (x[i] + 0.1)))
		{
			kx = i; break;
		}
	}
	for (j = 2; j < 4; j++)
	{
		if ((uy > (y[j] - 0.2)) && (uy <= (y[j] + 0.2)))
		{
			ry = j; break;
		}
	}
	if (tx <= (x[1] + 0.1)) kx = 1; if (tx > (x[4] - 0.1)) kx = 4;
	if (uy <= (y[1] + 0.2)) ry = 1; if (uy > (y[4] - 0.2)) ry = 4;

	p22 = 0.0;
	for (k = kx - 1; k <= kx + 1; k++)
	{
		temp1 = 1.0; temp3 = 0.0;
		for (t = kx - 1; t <= kx + 1; t++)
		{
			if (t == k) continue;
			temp1 *= (tx - x[t]) / (x[k] - x[t]);
		}

		for (r = ry - 1; r <= ry + 1; r++)
		{
			temp2 = 1;
			for (t = ry - 1; t <= ry + 1; t++)
			{
				if (t == r) continue;
				temp2 *= (uy - y[t]) / (y[r] - y[t]);
			}
			temp3 += temp2 * f[k][r];
		}
		p22 += temp1 * temp3;
	}
	return p22;
}

double power(double x, int y) //求X的y次方
{
	double temp = 1;
	for (int i = 1; i <= y; i++)
		temp *= x;
	return temp;
}


void approach(double x[], double y[], double fxy[][J], double C[I][I]) //二元拟合
{
	int i, j, r, k, s, t; //m=I=11,n=J=21,
	double error = 0.0;
	double B[I][I] = { 0.0 }, G[J][I] = { 0.0 }, P[I][J] = { 0.0 };
	double temp[I][I] = { 0.0 }, alfa[J][I] = { 0.0 }, gama[J][I] = { 0.0 }, b[I] = { 0.0 }; //此处为了方便alfa,gama均为转置形式

	for (i = 0; i < J; i++)
		for (j = 0; j < I; j++)
		{
			if (i < I) B[i][j] = 1.0;
			G[i][j] = 1.0;
		}
	outfile << "*******************************************" << endl;
	outfile << "选择过程中的k,delt值为：" << endl;
	for (k = 0;; k++) //M=N=k
	{
		error = 0.0;
		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
			{
				if (j < I) C[i][j] = 0.0;
				P[i][j] = 0.0;
			} //free error C P

		for (i = 0; i < I; i++)
			for (t = 1; t <= k; t++) B[i][k] *= x[i]; //构造φr,k=0,1,2...
		for (i = 0; i <= k; i++)
			for (j = 0; j <= k; j++)
				for (t = 0; t < I; t++) temp[i][j] += B[t][i] * B[t][j]; //B'*B

		for (j = 0; j < J; j++)
		{
			for (i = 0; i <= k; i++) 
				for (t = 0; t < I; t++) b[i] += B[t][i] * fxy[t][j];//B'*uj
			cout << endl << endl;
			Guass_choose(temp, b, alfa[j], k + 1); //alfa为转置
			for (i = 0; i <= k; i++) b[i] = 0.0;
		}
		for(auto j = 0;j < J;++j)
		{
			for(auto p = 0;p <= k;++p)
			{
				cout << alfa[j][p] << "\t";
			}
			cout << endl;
		}

		for (i = 0; i <= k; i++) //free temp
			for (t = 0; t <= k; t++) temp[i][t] = 0.0;

		for (i = 0; i < J; i++)
			for (t = 1; t <= k; t++) G[i][k] *= y[i]; //构造ψs,k=0,1,2...
		for (i = 0; i <= k; i++)
			for (j = 0; j <= k; j++)
				for (t = 0; t < J; t++) temp[i][j] += G[t][i] * G[t][j]; //G'*G

		for (j = 0; j < J; j++)
		{
			for (i = 0; i <= k; i++) b[i] = G[j][i]; //G'*ej
			Guass_choose(temp, b, gama[j], k + 1); //gama为转置
			for (i = 0; i <= k; i++) b[i] = 0.0;
		}

		for (i = 0; i <= k; i++) //free temp
			for (t = 0; t <= k; t++) temp[i][t] = 0.0;

		for (r = 0; r <= k; r++)
			for (s = 0; s <= k; s++)
				for (j = 0; j < J; j++) C[r][s] += alfa[j][r] * gama[j][s];

		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
			{
				for (r = 0; r <= k; r++)
					for (s = 0; s <= k; s++) P[i][j] += C[r][s] * power(x[i], r)*power(y[j], s);
				error += (P[i][j] - fxy[i][j])*(P[i][j] - fxy[i][j]);
			}
		outfile << k << "  " << error << endl;
		if (error < 1e-7) break;
	}
	outfile << "*******************************************" << endl;
	outfile << "达到精度的k,delt值为：" << endl;
	outfile << "k=" << k << "  " << "delt=" << error << endl;
	outfile << "*******************************************" << endl;
	outfile << "系数c[r][s]为：" << endl;
	for (r = 0; r <= k; r++)
	{
		outfile << "第" << r << "行:" << endl;
		for (s = 0; s <= k; s++)
		{
			outfile << C[r][s] << "  ";
			if ((s + 1) % 3 == 0) outfile << endl;
		}
	}
}

void main()
{
	double variable[L] = { 1,1,1,1 };//t,u,v,w
	double fz[6][6] = { 0.0 }, x[I] = { 0.0 }, y[J] = { 0.0 };
	double fxy[I][J] = { 0.0 }, C[I][I] = { 0.0 }, P[8][5] = { 0.0 };//f(x,y)
	int i, j, r, s, flag = 0;

	initial(fz);//z=f(u,t)
	outfile.open("output.txt");
	outfile << "*******************************************" << endl;
	outfile << "数表f(xi,yj)(i=0,1,...10;j=0,1,...20):" << endl;
	for (i = 0; i < I; i++)
	{
		for (j = 0; j < J; j++)
		{
			x[i] = 0.08*i; y[j] = 0.5 + 0.05*j;
			Newton(variable, x[i], y[j]);
			fxy[i][j] = Lagrange(variable[0], variable[1], fz);
			outfile << setprecision(0) << resetiosflags(ios::scientific);
			outfile << "f(" << x[i] << "," << y[j] << ")=";
			outfile << setprecision(12) << setiosflags(ios::scientific);
			outfile << fxy[i][j] << "      ";
			if (j % 2 == 1) outfile << endl;
		}
		outfile << endl << endl;
	}
	approach(x, y, fxy, C);

	outfile << "*******************************************" << endl;
	outfile << "数表f*(x,y):" << endl;
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 5; j++)
		{
			x[i] = 0.1*(i + 1); y[j] = 0.5 + 0.2*(j + 1);
			Newton(variable, x[i], y[j]);
			fxy[i][j] = Lagrange(variable[0], variable[1], fz);
			outfile << setprecision(0) << resetiosflags(ios::scientific);
			outfile << "f*(" << x[i] << "," << y[j] << ")=";
			outfile << setprecision(12) << setiosflags(ios::scientific);
			outfile << fxy[i][j] << "      ";
			if (j % 2 == 1) outfile << endl;
		}
		outfile << endl << endl;
	}
	outfile << "*******************************************" << endl;
	outfile << "数表p*(x,y):" << endl;
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 5; j++)
		{
			for (r = 0; r <= 5; r++)
				for (s = 0; s <= 5; s++)
					P[i][j] += C[r][s] * power(x[i], r)*power(y[j], s);
			outfile << setprecision(0) << resetiosflags(ios::scientific);
			outfile << "p*(" << x[i] << "," << y[j] << ")=";
			outfile << setprecision(12) << setiosflags(ios::scientific);
			outfile << P[i][j] << "      ";
			if (j % 2 == 1) outfile << endl;
		}
		outfile << endl << endl;
	}
}
