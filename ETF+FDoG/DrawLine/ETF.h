#pragma once

#include <math.h>
#include <cv.h>

class ETF
{
private:
	uchar* gray;
	double* gradient;
	double* gradientMag;
	double* tangent;
	double max_grad;
	int width;
	int height;
	int wideStep;
	int iterations;
	int radius;
public:
	~ETF();
	void Init(uchar* img, double* tan, int w, int h, int it, int r);
	void ETFTransfer();
private:
	void GetTangent();
	void RefineTangent();
	void RefineTangent2();
	inline double DistanceSquare(int x1, int y1, int x2, int y2);
	inline double Descartes(double x1, double y1, double x2, double y2);
	inline double Norm(double x, double y);
	inline void SobelFilter(double ul,  double um, double ur, double ml, double mr, double ll, double lm, double lr, double* x, double* y, double* mag, double* max_grad);
	inline void Normalize(double* x, double* y);
};