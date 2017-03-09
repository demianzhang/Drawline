#include "stdafx.h"
#include "ETF.h"
#include <fstream>
using namespace std;

#define ERROR_X 1024
#define ERROR_Y 1024

ETF::~ETF()
{
	delete[] gradient;
	delete[] gradientMag;
}
void ETF::Init(uchar* img, double* tan, int w, int h, int it, int r)
{
	gray = img;
	tangent = tan;
	width = w;
	height = h;
	wideStep = 2 * width;
	gradient = new double[wideStep * height];//Calculate gradient of the gray scale image
	gradientMag = new double[width * height];//gradient magnitude
	max_grad = -1.0;
	iterations= it;
	radius = r;
}

void ETF::ETFTransfer()
{
	GetTangent();
	RefineTangent();
}

void ETF::GetTangent()
{
	max_grad = -1.0;
	//计算图像每一点的梯度场和归一化的梯度值
	for(unsigned int j = 1; j < height - 1; j++)//get gradient of gradient map
	{
		for(unsigned int i = 1; i < width - 1; i++)
		{
			SobelFilter(gray[(j + 1)*width + i - 1], gray[(j + 1)*width + i], gray[(j + 1)*width + i + 1],
				gray[j*width + i - 1], gray[j*width + i + 1], gray[(j - 1)*width + i - 1], gray[(j - 1)*width + i], gray[(j - 1)*width + i + 1],
				&gradient[j*wideStep + i*2], &gradient[j*wideStep + i*2 + 1], &gradientMag[j*width + i], &max_grad);
			//gradient[j*wideStep + i*2], gradient[j*wideStep + i*2 + 1]存储该点横向及纵向的亮度差分
		}
	}
	//求出与该点梯度向量垂直的流场方向向量
	for(unsigned int i = 1; i < width - 1; i++)
	{
		gradient[i*2] = gradient[wideStep + i*2];//top row取紧邻的下一行值
		gradient[i*2 + 1] = gradient[wideStep + i*2 + 1];
		gradientMag[i] = gradientMag[width + i];
		gradient[(height - 1)*wideStep + i*2] = gradient[(height - 2)*wideStep + i*2];//bottom row取紧邻的上一行值
		gradient[(height - 1)*wideStep + i*2 + 1] = gradient[(height - 2)*wideStep + i*2 + 1];
		gradientMag[(height - 1)*width + i] = gradientMag[(height - 2)*width + i];
	}
	for(unsigned int j = 1; j < height - 1; j++)
	{
		gradient[j*wideStep] = gradient[j*wideStep + 2];//left most column取紧邻的右列值
		gradient[j*wideStep + 1] = gradient[j*wideStep + 3];
		gradientMag[j*width]= gradientMag[j*width + 1];
		gradient[j*wideStep + (width - 1)*2] = gradient[j*wideStep + (width -2)*2];//right most column取紧邻的左列值
		gradient[j*wideStep + (width - 1)*2 + 1] = gradient[j*wideStep + (width -2)*2 + 1];
		gradientMag[j*width + width - 1] = gradientMag[j*width + width - 2];
	}
	gradient[0] = (gradient[2] + gradient[wideStep]) / 2;//upper left corner(右+下)/2
	gradient[1] = (gradient[3] + gradient[wideStep + 1]) / 2;
	gradientMag[0] = (gradientMag[1] + gradientMag[width]) / 2;
	gradient[wideStep - 2] = (gradient[wideStep - 4] + gradient[2*wideStep - 2]) / 2;//upper right corner(左+下)/2
	gradient[wideStep - 1] = (gradient[wideStep - 3] + gradient[2*wideStep - 1]) / 2;
	gradientMag[width - 1] = (gradientMag[width - 2] + gradientMag[2*width - 1]) / 2;
	gradient[wideStep * (height - 1)] = (gradient[wideStep * (height - 2)] + gradient[wideStep * (height - 1) + 2]) / 2;//bottom left corner(上+右)/2
	gradient[wideStep * (height - 1) + 1] = (gradient[wideStep * (height - 2) + 1] + gradient[wideStep * (height - 1) + 3]) / 2;
	gradientMag[width * (height - 1)] = (gradientMag[width * (height - 2)] + gradientMag[width * (height - 1) + 1]) / 2;
	gradient[wideStep * height - 2] = (gradient[wideStep * height - 4] + gradient[wideStep * (height - 1) - 2]) / 2;//bottom right corner(上+左)/2
	gradient[wideStep * height -1] = (gradient[wideStep * height - 3] + gradient[wideStep * (height - 1) - 1]) / 2;
	gradientMag[width * height - 1] = (gradientMag[width * height - 2] + gradientMag[width * (height - 1) - 1]) / 2;


	for(unsigned int j = 0; j < height; j++)//Calculate origin tangent from gradient(梯度与流场方向向量数量积为0解得)
	{
		for(unsigned int i = 0; i < width; i++)
		{
			tangent[j*wideStep + i*2] = -1.0f * gradient[j*wideStep + i*2 + 1];
			tangent[j*wideStep + i*2 + 1] = gradient[j*wideStep + i*2];
			gradientMag[j*width + i] /= max_grad;
		} 
	}
}

 
//迭代地根据算法定义的核函数对ETF进行平滑
void ETF::RefineTangent()
{
	double* tangentNew = gradient;
	double cx, cy, cg, symbol;
	double wm, wd, weight, Upx, Upy, Down;
	int count = iterations;
	while(count > 0)
	{
		count--;
		for(int j = 0; j < height; j++)//horizontal
		{
			for(int i = 0; i < width; i++)
			{
				cx = tangent[j*wideStep + i*2];
				cy = tangent[j*wideStep + i*2 + 1];
				cg = gradientMag[j*width + i];
				Upx = Upy = 0.0;
				Down = 1e-10;

				for(int r = -radius; r <= radius; r++)
				{	
					int x = i+r;
					int y = j;
					if(x < 0)
						x = 0;
					if(x >= width)
						x = width - 1;
					//wm = (tanh(gradientMag[y*width + x] - cg) + 1) / 2;
					wm = ((gradientMag[y*width + x] - cg) + 1) / 2;//控制梯度值之差的影响,y与x的梯度值相差越大，这一权值也越高
					wd = Descartes(cx, cy, tangent[y*wideStep + x*2], tangent[y*wideStep + x*2 + 1]);
					symbol = (wd >= 0.0f) ? 1.0f : -1.0f;//空间距离权值
					wd *= symbol;//控制方向的影响,梯度向量越趋近于平行这一权值越大，越趋近于垂直则越小
					weight = wm * wd;
					Down += weight;
					Upx += (weight * symbol * tangent[y*wideStep + x*2]);//Upx=tangent(cur)*wm*wd*symbol
					Upy += (weight * symbol * tangent[y*wideStep + x*2 + 1]);//Upy=tangent(cur)*wm*wd*symbol
				}
				Upx /= Down;
				Upy /=Down;
				Normalize(&Upx, &Upy);//Normalize
				tangentNew[j*wideStep + i*2] = Upx;
				tangentNew[j*wideStep + i*2 + 1] = Upy;
			}
		}
		for(int i = 0; i < wideStep * height; i++)//copy from new tangent to current tangent
			tangent[i] = tangentNew[i];

		for(int j = 0; j < height; j++)//vertical
		{
			for(int i = 0; i < width; i++)
			{
				cx = tangent[j*wideStep + i*2];
				cy = tangent[j*wideStep + i*2 + 1];
				cg = gradientMag[j*width + i];
				Upx = Upy = 0.0;
				Down = 1e-10;

				for(int r = -radius; r <= radius; r++)
				{	
					int x = i;
					int y = j+r;
					if(y < 0)
						y = 0;
					if(y >= height)
						y = height - 1;
					//wm = (tanh(gradientMag[y*width + x] - cg) + 1) / 2;
					wm = ((gradientMag[y*width + x] - cg) + 1) / 2;
					wd = Descartes(cx, cy, tangent[y*wideStep + x*2], tangent[y*wideStep + x*2 + 1]);
					symbol = (wd >= 0.0f) ? 1.0f : -1.0f;
					wd *= symbol;
					weight = wm * wd;
					Down += weight;
					Upx += (weight * symbol * tangent[y*wideStep + x*2]);
					Upy += (weight * symbol * tangent[y*wideStep + x*2 + 1]);
				}
				Upx /= Down;
				Upy /=Down;
				Normalize(&Upx, &Upy);//Normalize
				tangentNew[j*wideStep + i*2] = Upx;
				tangentNew[j*wideStep + i*2 + 1] = Upy;
			}
		}
		for(int i = 0; i < wideStep * height; i++)//copy from new tangent to current tangent
			tangent[i] = tangentNew[i];
	}
}

void ETF::RefineTangent2()
{
	double* tangentNew = gradient;
	//int radius = 5;// to change
	double cx, cy, cg, symbol;
	double wm, wd, weight, Upx, Upy, Down;
	int count = iterations;
	while(count > 0)
	{
		count--;
		for(int j = 0; j < height; j++)
		{
			for(int i = 0; i < width; i++)
			{
				cx = tangent[j*wideStep + i*2];
				cy = tangent[j*wideStep + i*2 + 1];
				cg = gradientMag[j*width + i];
				Upx = Upy = 0.0;
				Down = 1e-10;

				for(int k = -radius; k <= radius; k++)
				{
					int y = j + k;
					if(y < 0)
						y = 0;
					if(y >= height)
						y = height - 1;
					for(int r = -radius; r <= radius; r++)
					{	
						int x = i + r;
						if(x < 0)
							x = 0;
						if(x >= width)
							x = width - 1;
						wm = (tanh(gradientMag[y*width + x] - cg) + 1) / 2;
						//wm = gradientMag[y*width + x] - cg;
						wd = Descartes(cx, cy, tangent[y*wideStep + x*2], tangent[y*wideStep + x*2 + 1]);
						symbol = (wd >= 0.0f) ? 1.0f : -1.0f;
						wd *= symbol;
						weight = wm * wd;
						//Down += weight;
						Upx += (weight * symbol * tangent[y*wideStep + x*2]);
						Upy += (weight * symbol * tangent[y*wideStep + x*2 + 1]);
					}
				}
				//Upx /= Down;
				//Upy /=Down;
				Normalize(&Upx, &Upy);//Normalize
				tangentNew[j*wideStep + i*2] = Upx;
				tangentNew[j*wideStep + i*2 + 1] = Upy;
			}
		}
		for(int i = 0; i < wideStep * height; i++)//copy from new tangent to current tangent
			tangent[i] = tangentNew[i];
	}
}

inline double ETF::DistanceSquare(int x1, int y1, int x2, int y2)
{
	double dx = x1 - x2;
	double dy = y1 - y2;
	return (dx * dx + dy * dy);
}

inline double ETF::Descartes(double x1, double y1, double x2, double y2)
{
	if(x1 == 0.0 && y1 == 0.0 && x2 == 0.0 && y2 == 0.0)
		return 1.0;
	return (x1 * x2 + y1 * y2);
}

inline double ETF::Norm(double x, double y)
{
	double norm;
	norm = x*x + y*y;
	if(norm < 1e-5)
		return 1e-5;
	else
		return sqrt(norm);
}

/*计算图像的梯度场和归一化的梯度值*/
//一阶导数法：梯度算子
inline void ETF::SobelFilter(double ul, // upper left
	double um, // upper middle
	double ur, // upper right
	double ml, // middle left
	double mr, // middle right
	double ll, // lower left
	double lm, // lower middle
	double lr, // lower right
	double* x,
	double* y,
	double* mag, double* max_grad)
{
	//算子包含两组3x3的矩阵，分别为横向及纵向，将之与图像作平面卷积，即可分别得出横向及纵向的亮度差分近似值
	double gx = (double)(ur + 2*mr + lr - ul - 2*ml - ll) / (4.0);// * 255
	double gy = (double)(ul + 2*um + ur - ll - 2*lm - lr) / (4.0);// * 255

	double norm = sqrt(gx * gx + gy * gy);
	*mag = norm;
	if(*max_grad < norm)// get the greatest gradient
		*max_grad = norm;
	if(norm < 1e-10)
	{
		*x = 0.0f;
		*y = 1.0f;
	}
	else
	{
		*x = gx / norm;
		*y = gy / norm;
	}
}

inline void ETF::Normalize(double* x, double* y)
{
	double nx = *x;
	double ny = *y;
	double norm = sqrt(nx*nx + ny*ny);

	if(norm < 1e-10)
	{
		*x = 0.0f;
		*y = 1.0f;
	}
	else
	{
		*x = nx / norm;
		*y = ny / norm;
	}
}

