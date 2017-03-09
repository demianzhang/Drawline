#pragma once
#include <math.h>
#include <cv.h>

#define PI 3.14
class FDoG
{
	
private:
	uchar* image;//�ο�ͼ��
	uchar* result;//���ͼ��
	double* tangent;//ETF
	int width;
	int height;
	int wideStep;
	double sigmaC;//�����߿�
	double sigmaS;//
	double sigmaM;//��������������
	double noise;
	double threshold;
	int deleteline;//����ͼƬ�򻯵ĳ̶�
	int extend;//�����沿���������������̶�
public:
	~FDoG();
	void Init(uchar* img, uchar* res, double* tan, int w, int h, double sc, double sm, double p, double t,int d,int e);//��ʼ������
	void FDoGTransfer();//����FDoG����
	void DecreaseLine();//���м򻯲���
	void DFS(int i,int j);
	void _DFS(int i,int j);
	void kDFS(int i,int j);
	void DFS_(int i,int j);
	void prim(int u0);
	void CountEdge();
private:
	double* MakeGaussians(int* n, double sigma);//����sigma����һά��˹
	inline double Gaussian(double x, double mean, double sigma);
	inline int round(double x);
};

