#include "mat.hpp"
#include <iostream>
using namespace std;
using namespace lzh;

int main()
{
	Mat::setPrintW(7);
	Mat::setPrintType(FIXED);
	Mat::setPrintSignificantDigits(3);
	mat_t data[] = {
		-1, -1, -1,
		 0,  0,  0,
		 1,  1,  1
	};
	Mat kern_y(data, 3, 3);
	cout << "kern_y = \n" << kern_y << endl;
	Mat kern_x = (Mat_(3, 3) <<
		-1, 0, 1,
		-1, 0, 1,
		-1, 0, 1);
	cout << "kern_x = \n" << kern_x << endl;
	Mat x = zeros(3, 3);
	cout << "zeros(3, 3) = \n" << x << endl;
	Mat xx = range(-3, 3);
	cout << "range(-3, 3) = \n" << xx << endl;
	Mat x3 = linspace(-3, 3, 10);
	cout << "linspace(-3, 3, 10) = \n" << x3 << endl;
	Mat x4(kern_x, kern_y, ROW);
	cout << "Mat x4(kern_x, kern_y, ROW) = \n" << x4 << endl;
	Mat x5(x4, x4, COL);
	cout << "Mat(x4, x4, COL) = \n" << endl;
	x5.show(); cout << endl;
	Mat x1 = mRand(0, 10, Size(3, 3), true);
	cout << "x1 = mRand(0, 10, Size(3, 3), true) = \n" << x1 << endl;
	Mat x2 = mRand(0, 10, Size(3, 3), true);
	cout << "x2 = mRand(0, 10, Size(3, 3), true) = \n" << x2 << endl;
	cout << "x1.Det() = " << x1.Det() << endl;
	cout << "x1.t() = \n" << x1.t() << endl;
	cout << "x1.adj() = \n" << x1.adj() << endl;
	cout << "x1.inv() = \n" << x1.inv() << endl;
	cout << "x1.reshape(1, -1) = " << endl;
	Mat x_r = x1;
	x_r.reshape(1, -1);
	cout << x_r << endl;
	cout << "x1 + x2 = \n" << x1 + x2 << endl;
	cout << "x1 - x2 = \n" << x1 - x2 << endl;
	cout << "x1 * x2 = \n" << x1 * x2 << endl;
	cout << "x1 / x2 = \n" << x1 / x2 << endl;
	cout << "dot(x1, x2) = \n" << Dot(x1, x2) << endl;
	cout << "exp(x1) = \n" << mExp(x1) << endl;
	cout << "log(x1) = \n" << mLog(x1) << endl;
	cout << "mul(x1, x2) = \n" << Mult(x1, x2) << endl;
	cout << "sum(x1) = " << x1.sum() << endl;
	cout << "abs(x1) = \n" << mAbs(x1) << endl;
	cout << "mSqrt(x1) = \n" << mSqrt(x1) << endl;
	cout << "mPow(x1, 1/2) = \n" << mPow(x1, 1 / 2.0) << endl;
	cout << "mPow(x1, 2) = \n" << mPow(x1, 2) << endl;
	cout << "adj(x1) = \n" << adj(x1) << endl;
	Mat vec = Random(3);
	cout << "vec = Random(3) = \n" << vec << endl;
	cout << "norm(vec) = " << mNorm(vec) << endl;
	cout << "norm(vec, 2) = " << mNorm(vec, 2) << endl;
	tools::pause();
	return 0;
}