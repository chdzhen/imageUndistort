#pragma once
#include<iostream>
using namespace std;


namespace xx
{
	double p00 = 3.788e-05;
	double p10 = 0.8197;
	double p01 = -0.0001723;
	double p20 = -0.0001939;
	double p11 = -0.0003418;
	double p02 = -0.000223;
	double p30 = -0.02762;
	double p21 = 0.0009897;
	double p12 = -0.1976;
	double p03 = 0.0005008;
	double p40 = 0.0001498;
	double p31 = 0.0004051;
	double p22 = 0.0005187;
	double p13 = 0.0004314;
	double p04 = 0.0001884;
	double p50 = 0.222;
	double p41 = -0.0008036;
	double p32 = 0.7629;
	double p23 = -0.001359;
	double p14 = 0.5242;
	double p05 = -0.0003189;



	/*double p00;
	double 	p10;
	double 	p01;
	double 	p20;
	double 	p11;
	double 	p02;
	double 	p30;
	double 	p21;
	double 	p12;
	double 	p03;
	double 	p40;
	double 	p31;
	double 	p22;
	double 	p13;
	double	p04;
	double	p50;
	double	p41;
	double	p32;
	double	p23;
	double	p14;
	double	p05;*/


	double get_answer(double x, double y)
	{
		double answer = p00 + p10*x + p01*y + p20*pow(x, 2) + p11*x*y + p02*pow(y, 2) + p30*pow(x, 3) + p21*pow(x, 2)  * y
			+ p12*x*pow(y, 2) + p03*pow(y, 3) + p40*pow(x, 4) + p31*pow(x, 3)  * y + p22*pow(x, 2)  * pow(y, 2)
			+ p13*x*pow(y, 3) + p04*pow(y, 4) + p50*pow(x, 5) + p41*pow(x, 4)  * y + p32*pow(x, 3)  * pow(y, 2)
			+ p23*pow(x, 2)  * pow(y, 3) + p14*x*pow(y, 4) + p05*pow(y, 5);
		return answer;
	}

}

namespace yy {

	double p00 = 3.733e-05;
	double p10 = -0.0001491;
	double p01 = 0.7512;
	double p20 = -0.0001886;
	double p11 = -0.0003298;
	double p02 = -0.0002247;
	double p30 = 0.0003741;
	double p21 = -0.162;
	double p12 = 0.001015;
	double p03 = -0.131;
	double p40 = 0.0001473;
	double p31 = 0.0003907;
	double p22 = 0.0004968;
	double p13 = 0.0004164;
	double p04 = 0.0001976;
	double p50 = -0.0002079;
	double p41 = 0.3097;
	double p32 = -0.001251;
	double p23 = 0.7631;
	double p14 = -0.0009004;
	double p05 = 0.3904;



	/*double	p00;
	double	p10;
	double	p01;
	double	p20;
	double	p11;
	double	p02;
	double	p30;
	double	p21;
	double	p12;
	double	p03;
	double	p40;
	double	p31;
	double	p22;
	double	p13;
	double	p04;
	double	p50;
	double	p41;
	double	p32;
	double	p23;
	double	p14;
	double	p05;*/

	double get_answer(double x, double y)
	{
		double answer = p00 + p10*x + p01*y + p20*pow(x, 2) + p11*x*y + p02*pow(y, 2) + p30*pow(x, 3) + p21*pow(x, 2)  * y
			+ p12*x*pow(y, 2) + p03*pow(y, 3) + p40*pow(x, 4) + p31*pow(x, 3)  * y + p22*pow(x, 2)  * pow(y, 2)
			+ p13*x*pow(y, 3) + p04*pow(y, 4) + p50*pow(x, 5) + p41*pow(x, 4)  * y + p32*pow(x, 3)  * pow(y, 2)
			+ p23*pow(x, 2)  * pow(y, 3) + p14*x*pow(y, 4) + p05*pow(y, 5);
		return answer;
	}

}