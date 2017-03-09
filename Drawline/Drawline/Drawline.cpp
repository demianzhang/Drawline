// Drawline.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <vector>
#include <queue>
const int MAXN = 200;//最大连通分量个数
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
}xy[MAXN][100000];//xy用于记录各标记连通线条的点集

/**
* @brief 对输入图像进行细化
* @param src为输入图像,用cvThreshold函数处理过的8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白
* @param maxIterations限制迭代次数，如果不进行限制，默认为-1，代表不限制迭代次数，直到获得最终结果
* @return 为对src细化后的输出图像,格式与src格式相同，元素中只有0与1,1代表有元素，0代表为空白
*/
cv::Mat thinImage(const cv::Mat & src, const int maxIterations = -1)
{
	assert(src.type() == CV_8UC1);
	cv::Mat dst;
	int width = src.cols;
	int height = src.rows;
	src.copyTo(dst);
	int count = 0;  //记录迭代次数
	while (true)
	{
		count++;
		if (maxIterations != -1 && count > maxIterations) //限制次数并且迭代次数到达
			break;
		std::vector<uchar *> mFlag; //用于标记需要删除的点
		//对点标记
		for (int i = 0; i < height; ++i)
		{
			uchar * p = dst.ptr<uchar>(i);
			for (int j = 0; j < width; ++j)
			{
				//如果满足四个条件，进行标记
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
						//标记
						mFlag.push_back(p + j);
					}
				}
			}
		}

		//将标记的点删除
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
		{
			**i = 0;
		}

		//直到没有点满足，算法结束
		if (mFlag.empty())
		{
			break;
		}
		else
		{
			mFlag.clear();//将mFlag清空
		}

		//对点标记
		for (int i = 0; i < height; ++i)
		{
			uchar * p = dst.ptr<uchar>(i);
			for (int j = 0; j < width; ++j)
			{
				//如果满足四个条件，进行标记
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
						//标记
						mFlag.push_back(p + j);
					}
				}
			}
		}

		//将标记的点删除
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
		{
			**i = 0;
		}

		//直到没有点满足，算法结束
		if (mFlag.empty())
		{
			break;
		}
		else
		{
			mFlag.clear();//将mFlag清空
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
			/////////////////////////////DFS遍历//////////////////////////////////////
			{ DFS_(i, j, dst), T++; cout << pos << endl; Tpos[T - 1] = pos; pos = 0; }//连通域标记
		}
	}
	
}



int main(int argc, char*argv[])
{
	//freopen("result.txt","w",stdout);
	//获取图像
	cv::Mat src = cv::imread("image.jpg", cv::IMREAD_GRAYSCALE);
	if (src.empty())
	{
		std::cout << "读取文件失败！" << std::endl;
		return -1;
	}
	cv::bitwise_not(src, src);
	//将原图像转换为二值图像
	cv::threshold(src, src, 128, 1, cv::THRESH_BINARY);
	//图片细化抽取骨架
	cv::Mat dst = thinImage(src);
	//显示图像
	dst = dst * 255;
	//二值图颜色反转
	cv::bitwise_not(dst, dst);
	//保存细化后的图片
	cv::imwrite("res.jpg", dst);
	//DFS顺序存储连通域内的点
	StorePoint(dst);
////////////////////////绘图操作/////////////////////////////////////
	/*cv::imshow("dst1", dst);
	cv::waitKey(0);*/

	IplImage* res_img;
	//设置绘图背景图片
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
	    //延时1ms动态刷新
		cvWaitKey(1);
		
	}
	cvSaveImage("result.jpg", res_img);
	cvWaitKey(0);
	cvDestroyWindow("Result of ETF+FDoG");
	cvReleaseImage( &res_img);
//////////////////////////////////////////////////////////////
	return 0;
}

