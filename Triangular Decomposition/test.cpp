#include "Triangular Decomposition.h"
#include <iostream>
using namespace std;

int main()
{
	Eigen::MatrixXd A,L,U;
	int n;
	cin >> n;
	A.resize(n, n);
	for (auto i = 0; i < n; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			cin >> A(i, j);
		}
	}
	DoolittleLU(A, L, U);
	cout << L << endl << U << endl;
	CroutLU(A, L, U);
	cout << L << endl << U << endl;
	system("pause");
}