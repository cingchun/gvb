#include <iostream>

double qingchun::fun(const double _q4, const double _q3, const double _q2, const double _q1, const double _theta)
{
	double sinx = std::sin(_theta), cosx = std::cos(_theta),
		sinx2 = sinx*sinx, sinx3 = sinx2*sinx, sinx4 = sinx3*sinx,
		cosx2 = cosx*cosx, cosx3 = cosx2*cosx;
	return _q4*sinx4 + _q3*sinx3*cosx + _q2*sinx2*cosx2 + _q1*sinx*cosx3;
}
double qingchun::dfun(const double _q4, const double _q3, const double _q2, const double _q1, const double _theta)
{
	double sinx = std::sin(_theta), cosx = std::cos(_theta),
		sinx2 = sinx*sinx, sinx3 = sinx2*sinx, sinx4 = sinx3*sinx,
		cosx2 = cosx*cosx, cosx3 = cosx2*cosx, cosx4 = cosx3*cosx;
	return  -_q3*sinx4 + (4.0*_q4 - 2.0*_q2)*sinx3*cosx + 3.0*(_q3 - _q1)*sinx2*cosx2 + 2.0*_q2*sinx*cosx3 + _q1*cosx4;
}
double qingchun::ddfun(const double _q4, const double _q3, const double _q2, const double _q1, const double _theta)
{
	double sinx = std::sin(_theta), cosx = std::cos(_theta),
		sinx2 = sinx*sinx, sinx3 = sinx2*sinx, sinx4 = sinx3*sinx,
		cosx2 = cosx*cosx, cosx3 = cosx2*cosx, cosx4 = cosx3*cosx;
	return (2.0*_q2 - 4.0*_q4)*sinx4 + (6.0*_q1 - 10.0*_q3)*sinx3*cosx + 12.0*(_q4 - _q2)*sinx2*cosx2 + (6.0*_q3 - 10.0*_q1)*sinx*cosx3 + 2.0*_q2*cosx4;
}
void qingchun::findroot(const double _q4, const double _q3, const double _q2, const double _q1, const double _q0, vector_1D_d &root, const double _precision)
{
	double  q4 = 1.0, q3 = _q3 / _q4, q2 = _q2 / _q4, q1 = _q1 / _q4, q0 = _q0 / _q4;
	vector_1D_d _root; // 存 y root
	findroot(-1.0, q2, 4.0*q0 - q3*q1, q0*q3*q3 + q1*q1 - 4.0*q0*q2, _root, _precision);// find y
	for (int i = 0; i < _root.size(); ++i)
	{
		double f = 0.5*q3*_root[i] - q1, d = std::sqrt(0.25*q3*q3 - q2 + _root[i]);
		findroot(1.0, q3 / 2.0 - d, 0.5*(_root[i] - f / d), root, _precision);
		findroot(1.0, q3 / 2.0 + d, 0.5*(_root[i] + f / d), root, _precision);
	}
}
void qingchun::findroot(const double _q3, const double _q2, const double _q1, const double _q0, vector_1D_d &root, const double _precision)
{
	// 方法一： 谢国芳
	// D = _q2*_q2 - 3.0*_q3*_q1;
	//if (std::fabs(D - 0.0) < precision) // D=0
	//{
	//	double c = _q2*_q2*_q2 - 27.0*_q3*_q3*_q0;
	//	if (std::fabs(c - 0.0) < precision) root.push_back(-_q2 / (3.0*_q3)); // 3 重根
	//	else root.push_back((-_q2 + std::pow(c, 1.0 / 3.0)) / (3.0*_q3));// 1 实根
	//}
	//else// D!=0
	//{
	//	double r = (9.0*_q3*_q2*_q1 - 2.0*_q2*_q2*_q2 - 27.0*_q3*_q3*_q0) / (2.0*std::sqrt(std::pow(std::fabs(D), 3.0)));// key ratio
	//	if (D<0.0) // D<0
	//	{
	//		double K = std::pow(r + std::sqrt(r*r + 1), 1.0 / 3.0);
	//		root.push_back((-_q2 + std::sqrt(std::fabs(D))*(K - 1.0 / K)) / (3.0*_q3));// 1 实根
	//	}
	//	else //D>0
	//	{
	//		double K = std::pow(r + std::sqrt(r*r - 1), 1.0 / 3.0);
	//		if (std::fabs(std::fabs(r) - 1.0)<precision)// |r|=1 2重实根和另一个实根
	//		{
	//			if (std::fabs(r - 1.0)<precision) // r=1  
	//			{
	//				root.push_back((-_q2 + 2.0*std::sqrt(D)) / (3.0*_q3));
	//				root.push_back((-_q2 - std::sqrt(D)) / (3.0*_q3));// 2 重实根
	//			}
	//			else // r=-1
	//			{
	//				root.push_back((-_q2 - 2.0*std::sqrt(D)) / (3.0*_q3));
	//				root.push_back((-_q2 + std::sqrt(D)) / (3.0*_q3));// 2 重实根
	//			}
	//		}
	//		else if (std::fabs(r) - 1.0>precision) // |r|>1 
	//		{
	//			root.push_back((-_q2 + std::sqrt(D)*(K + 1.0 / K)) / (3.0*_q3)); // 1 实根
	//		}
	//		else // |r|<1 3 互异实根
	//		{
	//			double theta = 1.0 / std::cos(r);
	//			root.push_back((-_q2 + 2.0*std::sqrt(D)*std::cos(theta / 3.0)) / (3.0*_q3));
	//			root.push_back((-_q2 + 2.0*std::sqrt(D)*std::cos((theta + 2.0 * PI) / 3.0)) / (3.0*_q3));
	//			root.push_back((-_q2 + 2.0*std::sqrt(D)*std::cos((theta - 2.0*PI) / 3.0)) / (3.0*_q3));
	//		}
	//	}
	//}

	// 方法二
	double _q32 = _q3*_q3, _q33 = _q32*_q3, _q22 = _q2*_q2, _q23 = _q22*_q2,
		p = (3.0*_q3*_q1 - _q22) / (3.0*_q32), q = (2.0*_q23 - 9.0*_q1*_q2*_q3 + 27.0*_q32*_q0) / (27.0*_q33), c = q*q / 4.0 + p*p*p / 27.0;
	if (std::fabs(c - 0.0)<_precision) // c=0
	{
		root.push_back(-2.0*std::sqrt(-p / 3.0) - _q2 / (3.0*_q3));
		root.push_back(std::sqrt(-p / 3.0) - _q2 / (3.0*_q3));
	}
	else if (c>0) // c>0
	{
		root.push_back(std::pow(-q / 2.0 + std::sqrt(c), 1.0 / 3.0) + std::pow(-q / 2.0 - std::sqrt(c), 1.0 / 3.0) - _q2 / (3.0*_q3));
	}
	else
	{
		double theta = std::acos((-q*std::sqrt(-27.0*p)) / (2.0*p*p));
		root.push_back(2.0*std::sqrt(-p / 3.0)*std::cos(theta / 3.0) - _q2 / (3.0*_q3));
		root.push_back(2.0*std::sqrt(-p / 3.0)*std::cos((theta - 2.0*PI) / 3.0) - _q2 / (3.0*_q3));
		root.push_back(2.0*std::sqrt(-p / 3.0)*std::cos((theta + 2.0*PI) / 3.0) - _q2 / (3.0*_q3));
	}
}
void qingchun::findroot(const double _q2, const double _q1, const double _q0, vector_1D_d &root, const double _precision)
{
	double c = _q1*_q1 - 4.0*_q0*_q2;
	if (std::fabs(c - 0.0) < _precision) root.push_back(-_q1 / (2.0*_q2)); // c=0
	else if (c > 0) // c>0
	{
		root.push_back((-_q1 + std::sqrt(c)) / (2.0*_q2));
		root.push_back((-_q1 - std::sqrt(c)) / (2.0*_q2));
		return;
	}
	else
	{
		//std::cerr << boost::format(" Warning: delta < 0 unsolvable\n");
		return;
	}
}
void qingchun::findroot(const double _q1, const double _q0, vector_1D_d &root, const double _precision)
{
	if (std::fabs(_q1-0.0)>=_precision)
	{
		root.push_back(-_q0 / _q1);
		return;
	}
	else
	{
		//std::cerr << boost::format(" Warning: _q1 < 0 unsolvable\n");
		return;
	}
}
double qingchun::findroot_theta(const double _q4, const double _q3, const double _q2, const double _q1, const double _precision)
{
	std::vector<double> _root;
	double q4 = _q3, q3 = 2.0*(_q2 - 2.0*_q4), q2 = 3.0*(_q1 - _q3), q1 = -2.0*_q2, q0 = -_q1;// f(x) x=tan(theta) 
	if (std::fabs(q4 - 0.0)>_precision) // quartic equation
	{
		findroot(q4, q3, q2, q1, q0, _root, _precision);
	}
	else if (std::fabs(q3 - 0.0)>_precision) // cubic equation
	{
		findroot(q3, q2, q1, q0, _root, _precision);
	}
	else if (std::fabs(q2 - 0.0)>_precision) // quadratic equation
	{
		findroot(q2, q1, q0, _root, _precision);
	}
	else if (std::fabs(q1 - 0.0)>_precision) // linear equation
	{
		findroot(q1, q0, _root, _precision);
	}
	else // constant equation 
	{
		return 0; // 无意义
	}
	// theta=arctan(x) 
	if (_root.size() == 0) // 无解时， 输出 0
	{
		return 0;
	}
	else// 有解时， 输出最小的那个根
	{
		double root = 0, f_min = 0, _f_min;
		for (int i = 0; i < _root.size(); ++i)
		{
			_root[i] = std::atan(_root[i]);
			_f_min = fun(_q4, _q3, _q2, _q1, _root[i]);
			if (_f_min<f_min)
			{
				root = _root[i];
				f_min = _f_min;
			}
		}
		return root;
	}
}
