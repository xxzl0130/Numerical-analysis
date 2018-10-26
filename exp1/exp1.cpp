#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
using namespace std;

typedef double Float;

#define R		2		//�°����
#define S		2		//�ϰ����
#define DIM		501		//ά��
#define eps		1e-12	//�����

Float A[R + S + 1][DIM], B[R + S + 1][DIM];

/**
 * \brief ��ʼ������
 */
void init();
Float powerMethod();
Float inversePowerMethod();
void normalize(Float *vec, size_t n);
void LU();

int main()
{
	init();

	memcpy_s(B, sizeof B, A, sizeof A);
	const auto lambda1 = powerMethod();//�ݷ���ģ��������ֵ
	system("pause");
	memcpy_s(B, sizeof B, A, sizeof A);
	for (auto i = 0; i < DIM; ++i)//��ֵ����B=��A-lamk*I��
	{
		B[S][i] -= lambda1;
	}
	const auto lambda2 = powerMethod() + lambda1;//��ԭ��λ�Ƶ��ݷ���λ�ƺ�ģ��������ֵ
	const auto maxLambda = max(lambda1, lambda2);
	const auto minLambda = min(lambda1, lambda2);
	printf("lam1   = %12.12e\n", minLambda);
	printf("lam501 = %12.12e\n", maxLambda);

	memcpy_s(B, sizeof B, A, sizeof A);
	const auto lam3 = inversePowerMethod();
	printf("lams   = %12.12e\n", lam3);

	for (auto k = 1; k <= 39; ++k)//ѭ������u�ӽ�������ֵ
	{
		auto u = minLambda + k * (maxLambda - minLambda) / 40;
		//��ֵB=(A-miu*I)
		memcpy_s(B, sizeof B, A, sizeof A);
		for(auto i = 0;i < DIM;++i)
		{
			B[2][i] -= u;
		}
		const auto lambda = inversePowerMethod() + u;//�ô�ԭ��ķ��ݷ�����miu�ӽ�������ֵ
		printf("lami%02d = %12.12e\n", k, lambda);
	}

	const auto cond2 = abs(maxLambda / lam3);
	printf("cond2  = %12.12e\n", cond2);

	memcpy_s(B, sizeof B, A, sizeof A);
	LU();
	Float det = 1;
	for(auto i = 0;i < DIM;++i)
	{
		det *= B[2][i];
	}
	printf("det(A) = %12.12e\n", det);

	system("pause");
}

void init()
{
	for (auto i = 0; i < DIM; ++i)
	{
		const auto j = i + 1;
		A[2][i] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j);
		A[3][i] = 0.16;
		A[4][i] = -0.064;
	}
	for (auto i = 2; i < DIM; ++i)
	{
		A[0][i] = -0.064;
	}
	for (auto i = 1; i < DIM; ++i)
	{
		A[1][i] = 0.16;
	}
	A[3][500] = 0;
	A[4][500] = 0;
	A[4][499] = 0;
}

/**
 * \brief �ݷ�
 * \return ��ģ�������ֵ
 */
Float powerMethod()
{
	Float u[DIM] = {};
	Float y[DIM] = {};
	Float beta1 = 0, beta2 = 0;
	Float err;
	int count = 0;
	//u������ʼ��
	for (auto& i : u)
	{
		i = 1;
	}

	do
	{
		normalize(u, DIM);
		memcpy_s(y, sizeof y, u, sizeof u);//����y=u
		memset(u, 0, sizeof u);//���u
		for (auto i = 2; i <= DIM + 1; ++i)//����u(k)
		{
			int j;
			if (i < 4) 
				j = 0;
			else 
				j = i - 4;
			for (; j <= min(i,DIM - 1); j++)
			{
				u[i - 2] += B[i - j][j] * y[j];
			}
		}

		for(auto i = 0;i < DIM;++i)
		{
			beta2 += y[i] * u[i];
		}
		err = abs(beta2 - beta1) / abs(beta2);
		beta1 = beta2;
		beta2 = 0;
		++count;
	} while (err > eps);
	printf("%d\n", count);
	return beta1;
}

/**
 * \brief ���ݷ�
 * \return ��ģ��С����ֵ
 */
Float inversePowerMethod()
{
	Float u[DIM];
	Float y[DIM] = {};
	Float d[DIM] = {};
	Float beta1 = 0, beta2 = 0;
	Float err;

	//u������ʼ��
	for (auto& i : u)
	{
		i = 1;
	}

	LU();

	do
	{
		normalize(u, DIM);
		memcpy_s(y, sizeof y, u, sizeof u);//����y=u
		memset(u, 0, sizeof u);//���u

		for(auto i = 0;i < DIM;++i)
		{
			Float sum = 0;
			for(auto j = max(1,i - 1);j <= i;++j)
			{
				sum += B[i - j + 3][j - 1] * d[j - 1];
			}
			d[i] = y[i] - sum;
		}

		for (auto i = DIM - 1; i >= 0; --i)
		{
			Float sum = 0;
			for (auto j = i + 1; j < min(i + 3, DIM); ++j)
			{
				sum += B[i - j + 2][j] * u[j];
			}
			u[i] = (d[i] - sum) / B[2][i];
		}

		for (auto i = 0; i < DIM; ++i)
		{
			beta2 += y[i] * u[i];
		}
		err = abs(beta2 - beta1) / abs(beta2);
		beta1 = beta2;
		beta2 = 0;
	} while (err > eps);
	
	return 1.0 / beta1;
}

/**
 * \brief ������һ��
 * \param vec ����ָ��
 * \param n ��������
 */
void normalize(Float* vec, size_t n)
{
	Float s = 0;
	for(auto i = 0u;i < n;++i)
	{
		s += vec[i] * vec[i];
	}
	s = sqrt(s);
	for (auto i = 0u; i < n; ++i)
	{
		vec[i] /= s;
	}
}

/**
 * \brief ��B��LU�ֽ�
 */
void LU()
{
	for (auto i = 3; i <= 4; ++i)//����k=1ʱ��Ӧ�ľ���Ԫ���ȸ�ֵ
	{
		B[i][0] = B[i][0] / B[2][0];
	}
	for (auto k = 2; k <= DIM; ++k)//��k=2��ʼ��U�е�Ԫ�ش浽ϵ��������
	{
		for (auto j = k; j <= min(DIM,k + 1); ++j)
		{
			Float sum = 0;
			for (auto t = max(max(k - 2,j - 2),1); t < k; ++t)
			{
				sum += B[k - t + 2][t - 1] * B[t - j + 2][j - 1];
			}
			B[k - j + 2][j - 1] = B[k - j + 2][j - 1] - sum;
		}
		if (k < DIM)
		{
			for (auto i = k + 1; i <= min(DIM,k + 2); i++)
			{
				Float sum = 0;
				for (auto j = max(max(i - 2,k - 2),1); j < k; ++j)
				{
					sum += B[i - j + 2][j - 1] * B[j - k + 2][k - 1];
				}
				B[i - k + 2][k - 1] = B[i - k + 2][k - 1] - sum;
				B[i - k + 2][k - 1] = B[i - k + 2][k - 1] / B[2][k - 1];
			}
		}
	}
}

