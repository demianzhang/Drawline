#pragma once

#include "cv.h"
#include "highgui.h"
#include "ETF.h"
#include "FDoG.h"

struct Parameters
{
	int smooth_times;
	int iteration_times;
	int radius;
	double sigmaC;//�����߿�
	double sigmaM;//��������������
	double noise;
	double threshold;
	int deleteline;//����ͼƬ�򻯵ĳ̶�
	int extend;//�����沿���������������̶�
};

class DrawLine
{
private:
	IplImage* src_img;
	IplImage* res_img;
	unsigned int width;
	unsigned int height;
	double* tangent;
	struct Parameters param;
public:
	DrawLine();
	~DrawLine();
	void SetParameters(int st, int it, int rr, double sc, double sm, double pp, double tt,int d,int e);
	void SetSrcImg(IplImage* img);
	IplImage* GetSrcImg();
	double* GetETFImg();
	IplImage* GetResultImg();
	void TransferImage();
};