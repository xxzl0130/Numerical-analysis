#include "Gaussian Elimination.h"
#include <iostream>
using namespace std;

int main()
{
	Eigen::MatrixXd A;
	int n;
	cin >> n;
	A.resize(n, n + 1);
	for(auto i = 0;i < n;++i)
	{
		for(auto j = 0;j <= n;++j)
		{
			cin >> A(i, j);
		}
	}
	auto X1 = SequenceGaussian(A);
	auto X2 = MainElementGaussian(A);
	cout << X1.transpose() << endl << X2.transpose() << endl;
	system("pause");
}