#include "stdafx.h"
#include "DrawLine.h"
#include <fstream>
using namespace std;

DrawLine::DrawLine()
{
	tangent = 0;
	src_img = 0;
	res_img = 0;
}
DrawLine::~DrawLine()
{
	if(tangent)
		delete[] tangent;
	if(res_img)
		cvReleaseImage(&res_img);
	if(src_img)
		cvReleaseImage(&src_img);
}

void DrawLine::SetParameters(int st, int it, int rr, double sc, double sm, double pp, double tt,int d,int e)
{
	param.smooth_times = st;
	param.iteration_times = it;
	param.radius = rr;
	param.sigmaC = sc;
	param.sigmaM = sm;
	param.noise = pp;
	param.threshold = tt;
	param.deleteline=d;
	param.extend=e;
}

void DrawLine::SetSrcImg(IplImage* img)
{
	src_img = img;
	width = img->width;
	height = img->height;
}

IplImage* DrawLine::GetSrcImg()
{
	return src_img;
}

double* DrawLine::GetETFImg()
{
	return tangent;
}

IplImage* DrawLine::GetResultImg()
{
	return res_img;
}

void DrawLine::TransferImage()//do transfer
{
	IplImage* gray_img = cvCreateImage(cvGetSize(src_img), 8, 1);//Convert BGR image to gray scale
	cvCvtColor(src_img, gray_img, CV_BGR2GRAY);
	uchar* gray = (uchar*)gray_img->imageData;



	int count = param.smooth_times;
	while(count > 0)
	{
		count--;
		cvSmooth(gray_img, gray_img, CV_GAUSSIAN, 3, 3);//gaussian smooth
	}

	tangent = new double[2 * width * height];//Calculate origin tangent from gradient
	ETF etf;
	etf.Init(gray, tangent, width, height, param.iteration_times, param.radius);
	etf.ETFTransfer();



	res_img = cvCreateImage(cvGetSize(src_img), 8, 3);
	FDoG fdog;
	fdog.Init(gray, (uchar*)res_img->imageData, tangent, width, height, param.sigmaC, param.sigmaM, param.noise, param.threshold,param.deleteline,param.extend);
	fdog.FDoGTransfer();
	fdog.DecreaseLine();//��
    fdog.CountEdge();//��ͨ���ǣ���¼����
	fdog.prim(1);

	//cvReleaseImage(&gray_img);
}
/*
ע������һֱû���ҵ�bug��
��fdog.CountEdge();fdog.prim(1);
���ж�ε���Ҳû����ȫ��ͨ��Ч�������ر��
�����ϻ���ȫ��ͨ�ģ�Audrey Hepburn.jpg����ͨЧ����lena.jpgҪ�úܶ�
*/