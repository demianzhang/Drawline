// main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "_FDoG.h"
#include "_DrawLine.h"
#include "_ETF.h"
#include <string.h>


int main()
{
	freopen("key.txt","r",stdin);
	
	DrawLine draw;
	IplImage* src_img;
	IplImage* result_img1;
	src_img = cvLoadImage("lena.jpg",1);
	draw.SetSrcImg(src_img);
	cvNamedWindow("Source Image", 1); 
	cvShowImage("Source Image", draw.GetSrcImg());
	cvWaitKey(0);
	cvDestroyWindow("Source Image");
	int x,y,i=0;
	//载入key.txt特征点
	while(~scanf("%d %d",&x,&y))
	{
		key[i][0]=x;
		key[i][1]=y;
		i++;
	}
	len=i;
////////////////////////测试点的分布/////////////////////////////
//	for(int i=0;i<len;i++)
//{
//CvPoint centerpoint;
//centerpoint.x=int(key[i][0]);
//centerpoint.y=int(key[i][1]);
//cvCircle( src_img, centerpoint ,3 , CV_RGB(0,255,0),1, 8, 3 );
//} 
//cvNamedWindow("image_test",CV_WINDOW_AUTOSIZE);
//cvSaveImage("image.jpg",src_img);
//cvShowImage("image_test",src_img);
//cvWaitKey(0);

//////////////////////////////////////////////////////////////////
	draw.SetParameters(1, 3, 5, 1.00, 3.00, 0.99, 0.50,300,30);
	draw.TransferImage();
	result_img1 = draw.GetResultImg();
	cvNamedWindow("Result of ETF+FDoG", 1);
	cvShowImage("Result of ETF+FDoG", result_img1);
	cvWaitKey(0);
    cvDestroyWindow("Result of Algorithm1");
	
	return 0;
}