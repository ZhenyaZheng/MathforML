#include "mymath.h"
#include <iostream>
using namespace std;
int main()
{
	vector<double> a(20), b(20);
	vector<int> ia(20), ib(20);
	for (int i = 0; i < 20; ++i)
	{
		a[i] = 0.3 * (i + 1);
		b[i] = 1.1 * (i + 1);
		ia[i] = i % 5;
		ib[i] = (34 - i) % 5;
	}
	auto dval = MyMath::MyMath::getTtest(a, b);
	auto ival = MyMath::MyMath::getChitest(MyMath::MyMath::generateDiscreteAttributesCategoryIntersection(ia, ib, 5, 5));
	cout << "dval = " << dval << endl;
	cout << "ival = " << ival << endl;
	return 0;
}