#include <cstdio>
#include <cmath>
#include <cstring>

double A[5][501], B[5][501];

/**
 * \brief ³õÊ¼»¯¾ØÕó
 */
void init();

int main()
{
	init();
}

void init()
{
	for(auto i = 0;i < 501;++i)
	{
		const auto j = i + 1;
		A[2][i] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j);
	}
	for(auto i = 2;i < 501;++i)
	{
		A[0][i] = A[4][i] = -0.064;
	}
	for (auto i = 1; i < 501; ++i)
	{
		A[1][i] = A[3][i] = 0.16;
	}
	memcpy_s(B, sizeof B, A, sizeof A);
}