/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef COLORENGINE_H
#define COLORENGINE_H

#include <vector>

class ColorEngine{
public:
	ColorEngine(){}
	~ColorEngine(){}

	struct color{
		float r;
		float g;
		float b;
	};
	struct HSL{
		double X;
		double Y;
		double Z;
	};	

public:
	static void colorHSV(unsigned n, std::vector<color> &colors);
	static void colorHeat(unsigned n, std::vector<color> &colors);
	static void HslToRgb(std::vector<double>& mag, std::vector<color> &colors);
	static double HueToRgb(double p, double q, double t);
};


#endif
