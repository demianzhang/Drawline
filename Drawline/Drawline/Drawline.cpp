// Drawline.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <vector>
#include <queue>
const int MAXN = 200;//�����ͨ��������
const int INF = 1000000;
using namespace std;

static int T = 1, visit[1000][1000], pos = 0, minint = 0x7fffffff;
static int Tpos[MAXN];
static  struct p
{
	int x;
	int y;
	p(){ x = 0, y = 0; }
	p(int i, int j):x(i), y(j){}
}xy[MAXN][100000];//xy���ڼ�¼�������ͨ�����ĵ㼯

/**
* @brief ������ͼ�����ϸ��
* @param srcΪ����ͼ��,��cvThreshold�����������8λ�Ҷ�ͼ���ʽ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
* @param maxIterations���Ƶ���������������������ƣ�Ĭ��Ϊ-1���������Ƶ���������ֱ��������ս��
* @return Ϊ��srcϸ��������ͼ��,��ʽ��src��ʽ��ͬ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
*/
cv::Mat thinImage(const cv::Mat & src, const int maxIterations = -1)
{
	assert(src.type() == CV_8UC1);
	cv::Mat dst;
	int width = src.cols;
	int height = src.rows;
	src.copyTo(dst);
	int count = 0;  //��¼��������
	while (true)
	{
		count++;
		if (maxIterations != -1 && count > maxIterations) //���ƴ������ҵ�����������
			break;
		std::vector<uchar *> mFlag; //���ڱ����Ҫɾ���ĵ�
		//�Ե���
		for (int i = 0; i < height; ++i)
		{
			uchar * p = dst.ptr<uchar>(i);
			for (int j = 0; j < width; ++j)
			{
				//��������ĸ����������б��
				//  p9 p2 p3
				//  p8 p1 p4
				//  p7 p6 p5
				uchar p1 = p[j];
				if (p1 != 1) continue;
				uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
				uchar p8 = (j == 0) ? 0 : *(p + j - 1);
				uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);
				uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);
				uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);
				uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);
				uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);
				uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);
				if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
				{
					int ap = 0;
					if (p2 == 0 && p3 == 1) ++ap;
					if (p3 == 0 && p4 == 1) ++ap;
					if (p4 == 0 && p5 == 1) ++ap;
					if (p5 == 0 && p6 == 1) ++ap;
					if (p6 == 0 && p7 == 1) ++ap;
					if (p7 == 0 && p8 == 1) ++ap;
					if (p8 == 0 && p9 == 1) ++ap;
					if (p9 == 0 && p2 == 1) ++ap;

					if (ap == 1 && p2 * p4 * p6 == 0 && p4 * p6 * p8 == 0)
					{
						//���
						mFlag.push_back(p + j);
					}
				}
			}
		}

		//����ǵĵ�ɾ��
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
		{
			**i = 0;
		}

		//ֱ��û�е����㣬�㷨����
		if (mFlag.empty())
		{
			break;
		}
		else
		{
			mFlag.clear();//��mFlag���
		}

		//�Ե���
		for (int i = 0; i < height; ++i)
		{
			uchar * p = dst.ptr<uchar>(i);
			for (int j = 0; j < width; ++j)
			{
				//��������ĸ����������б��
				//  p9 p2 p3
				//  p8 p1 p4
				//  p7 p6 p5
				uchar p1 = p[j];
				if (p1 != 1) continue;
				uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
				uchar p8 = (j == 0) ? 0 : *(p + j - 1);
				uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);
				uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);
				uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);
				uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);
				uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);
				uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);

				if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
				{
					int ap = 0;
					if (p2 == 0 && p3 == 1) ++ap;
					if (p3 == 0 && p4 == 1) ++ap;
					if (p4 == 0 && p5 == 1) ++ap;
					if (p5 == 0 && p6 == 1) ++ap;
					if (p6 == 0 && p7 == 1) ++ap;
					if (p7 == 0 && p8 == 1) ++ap;
					if (p8 == 0 && p9 == 1) ++ap;
					if (p9 == 0 && p2 == 1) ++ap;

					if (ap == 1 && p2 * p4 * p8 == 0 && p2 * p6 * p8 == 0)
					{
						//���
						mFlag.push_back(p + j);
					}
				}
			}
		}

		//����ǵĵ�ɾ��
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
		{
			**i = 0;
		}

		//ֱ��û�е����㣬�㷨����
		if (mFlag.empty())
		{
			break;
		}
		else
		{
			mFlag.clear();//��mFlag���
		}
	}
	return dst;
}


void DFS_(int i, int j, cv::Mat& dst)
{
	if (visit[i][j] != 0 || i < 0 || i >= dst.rows || j < 0 || j >= dst.cols)return;
	if ((unsigned)dst.at<uchar>(i, j) == 0)
	{
		visit[i][j] = T;
		xy[T][pos].x = i;
		xy[T][pos].y = j;
		pos++;
		DFS_(i - 1, j - 1, dst);
		DFS_(i - 1, j, dst);
		DFS_(i - 1, j + 1, dst);
		DFS_(i, j - 1, dst);
		DFS_(i, j + 1, dst);
		DFS_(i + 1, j - 1, dst);
		DFS_(i + 1, j, dst);
		DFS_(i + 1, j + 1, dst);
	}
}

queue<p> q;
void StorePoint(cv::Mat& dst)
{
	memset(visit, 0, sizeof(visit));
	minint = 0x7fffffff;
	T = 1;
	pos = 0;
	for (int i = 1; i < dst.rows - 1; i++)
	{
		for (int j = 1; j < dst.cols - 1; j++)
		{
			
			if ((unsigned)dst.at<uchar>(i, j) == 255 || visit[i][j] != 0)continue;
			if ((unsigned)dst.at<uchar>(i, j)==0&& visit[i][j] == 0)
			/////////////////////////////DFS����//////////////////////////////////////
			{ DFS_(i, j, dst), T++; cout << pos << endl; Tpos[T - 1] = pos; pos = 0; }//��ͨ����
		}
	}
	
}



int main(int argc, char*argv[])
{
	//freopen("result.txt","w",stdout);
	//��ȡͼ��
	cv::Mat src = cv::imread("image.jpg", cv::IMREAD_GRAYSCALE);
	if (src.empty())
	{
		std::cout << "��ȡ�ļ�ʧ�ܣ�" << std::endl;
		return -1;
	}
	cv::bitwise_not(src, src);
	//��ԭͼ��ת��Ϊ��ֵͼ��
	cv::threshold(src, src, 128, 1, cv::THRESH_BINARY);
	//ͼƬϸ����ȡ�Ǽ�
	cv::Mat dst = thinImage(src);
	//��ʾͼ��
	dst = dst * 255;
	//��ֵͼ��ɫ��ת
	cv::bitwise_not(dst, dst);
	//����ϸ�����ͼƬ
	cv::imwrite("res.jpg", dst);
	//DFS˳��洢��ͨ���ڵĵ�
	StorePoint(dst);
////////////////////////��ͼ����/////////////////////////////////////
	/*cv::imshow("dst1", dst);
	cv::waitKey(0);*/

	IplImage* res_img;
	//���û�ͼ����ͼƬ
	res_img = cvLoadImage("background.jpg", 1);

	cvNamedWindow("Result of ETF+FDoG", 1);
	for (int i = 1; i < T; i++)
	for (int j = 0; j < Tpos[i]; j++)
	{
		CvPoint centerpoint;
		centerpoint.x = xy[i][j].y;
		centerpoint.y = xy[i][j].x;
		cvCircle(res_img, centerpoint, 2, CV_RGB(0, 0, 0), -1, 8, 0);
		cvShowImage("Result of ETF+FDoG", res_img);
	    //��ʱ1ms��̬ˢ��
		cvWaitKey(1);
		
	}
	cvSaveImage("result.jpg", res_img);
	cvWaitKey(0);
	cvDestroyWindow("Result of ETF+FDoG");
	cvReleaseImage( &res_img);
//////////////////////////////////////////////////////////////
	return 0;
}

