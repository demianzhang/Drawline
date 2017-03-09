#pragma once
#include <math.h>
#include <cv.h>

#define PI 3.14
class FDoG
{
	
private:
	uchar* image;//参考图像
	uchar* result;//结果图像
	double* tangent;//ETF
	int width;
	int height;
	int wideStep;
	double sigmaC;//控制线宽
	double sigmaS;//
	double sigmaM;//控制线条连续性
	double noise;
	double threshold;
	int deleteline;//控制图片简化的程度
	int extend;//控制面部特征轮廓的清晰程度
public:
	~FDoG();
	void Init(uchar* img, uchar* res, double* tan, int w, int h, double sc, double sm, double p, double t,int d,int e);//初始化参数
	void FDoGTransfer();//进行FDoG操作
	void DecreaseLine();//进行简化操作
	void DFS(int i,int j);
	void _DFS(int i,int j);
	void kDFS(int i,int j);
	void DFS_(int i,int j);
	void prim(int u0);
	void CountEdge();
private:
	double* MakeGaussians(int* n, double sigma);//根据sigma生成一维高斯
	inline double Gaussian(double x, double mean, double sigma);
	inline int round(double x);
};

