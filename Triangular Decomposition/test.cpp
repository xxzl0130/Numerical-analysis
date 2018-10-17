#include "Triangular Decomposition.h"
#include <iostream>
using namespace std;

int main()
{
	Eigen::MatrixXd A,L,U,B;
	int n;
	cin >> n;
	A.resize(n, n);
	B.resize(n, 1);
	for (auto i = 0; i < n; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			cin >> A(i, j);
		}
		cin >> B(i);
	}
	DoolittleLU(A, L, U);
	cout << L << endl << U << endl;
	cout << "------------------" << endl;
	MainElementDoolittleLU(A, B, L, U);
	cout << L << endl << U << endl;
	cout << B << endl;
	system("pause");
}