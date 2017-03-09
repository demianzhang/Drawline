#include "stdafx.h"
#include "FDoG.h"
#include <fstream>
#include <iostream>
const int MAXN=200;//最大连通分量个数
const int INF=1000000;
using namespace std;

FDoG::~FDoG()
{

}
void FDoG::Init(uchar* img, uchar* res, double* tan, int w, int h, double sc, double sm, double p, double t,int d,int e)
{
	image = img;
	result = res;
	tangent = tan;
	width = w;
	height = h;
	wideStep = 2 * width;
	sigmaC = sc;
	sigmaS = 1.6 * sigmaC;
	sigmaM = sm;
	noise = p;
	threshold = t;
	deleteline=d;
	extend=e;
}

void FDoG::FDoGTransfer()
{
	int NC, NS, NM;//卷积核的大小
	double* gaussC = MakeGaussians(&NC, sigmaC);// get gaussians from sigma
	double* gaussS = MakeGaussians(&NS, sigmaS);
	double* gaussM = MakeGaussians(&NM, sigmaM);

	double* F = new double[width * height];
	int wc = NC - 1;
	int ws = NS - 1;
	int x, y;
	double d_x, d_y;
	double gx, gy;
	double sumC, sumS, weightC, weightS;
	
	for(int j = 0; j < height; j++)//Get integral along gradient（梯度）
	{
		for(int i = 0; i < width; i++)
		{
			sumC = sumS = 0.0;
			weightC = weightS = 1e-10;
			gy = tangent[j*wideStep + 2*i];
			gx = -1.0 * tangent[j*wideStep + 2*i + 1];
			if(gx == 0.0 && gy == 0.0)//handle (0,0) gradient
			{
				F[j*width + i] = 255.0 * (1.0 - noise);
				continue;
			}
			
			for(int k = -ws; k <= ws; k++)
			{
				d_x = (double)i + (double)k*gx;
				d_y = (double)j + (double)k*gy;
				x = round(d_x);
				y = round(d_y);
		
				if(x < 0)
					x = 0;
				if(x >= width)
					x = width - 1;
				if(y < 0)
					y = 0;
				if(y >= height)
					y = height - 1;
				if(abs(k) <= wc)
				{
					sumC += ((double)image[y*width + x] * gaussC[abs(k)]);
					weightC += gaussC[abs(k)];
				}
				sumS += ((double)image[y*width + x] * gaussS[abs(k)]);
				weightS += gaussS[abs(k)];
			}
			F[j*width + i] = (sumC / weightC) - noise * (sumS / weightS);
		//在每一点处沿着该点的梯度方向作一维DoG，并将结果保存
		}
	}
	//曲线S代表的是方向场的流向，为了得到更连续的线条，需要再沿着S进行一次高斯卷积
	int wm = NM - 1;
	double tx, ty;
	double sumM, weightM;
	double h;
	double stepSize = 1.0;
	
	for(int j = 0; j < height; j++)//Get integral along flow
	{
		for(int i = 0; i < width; i++)
		{
			sumM = F[j*width + i] * gaussM[0];
			weightM = 1e-10 + gaussM[0];
			x = i;
			y = j;
			d_x = (double)x;
			d_y = (double)y;
			for(int k = 1; k <= wm; k++)//positive direction tangent
			{
				tx = tangent[y*wideStep + 2*x];
				ty = tangent[y*wideStep + 2*x + 1];
				if(tx == 0.0 && ty == 0.0)//omit later point when encounters (0,0) tangent
					break;
				d_x += (tx * stepSize);//stepSize = 1
				d_y += (ty * stepSize);
				if(d_x < 0.0 || d_x >= width || d_y < 0.0 || d_y >= height)
					break;
				x = round(d_x);
				y = round(d_y);
				if(x < 0 || x >= width || y < 0 || y >= height)
					break;
				weightM += gaussM[k];
				sumM += (F[y*width + x] * gaussM[k]);
			}
			x = i;
			y = j;
			d_x = (double)x;
			d_y = (double)y;
			for(int k = 1; k <= wm; k++)//negative direction tangent
			{	
				tx = -1.0 * tangent[y*wideStep + 2*x];
				ty = -1.0 * tangent[y*wideStep + 2*x + 1];
				if(tx == 0.0 && ty == 0.0)//omit later point when encounters (0,0) tangent
					break;
				d_x += (tx * stepSize);//stepSize = 1
				d_y += (ty * stepSize);
				if(d_x < 0.0 || d_x >= width || d_y < 0.0 || d_y >= height)
					break;
				x = round(d_x);
				y = round(d_y);
				if(x < 0 || x >= width || y < 0 || y >= height)
					break;
				weightM += gaussM[k];
				sumM += F[y*width + x] * gaussM[k];
			}
			h = sumM / weightM;
		//根据设定的阈值 对图像进行二值化RGB
			if(h < 0 && (1.0 + tanh(h)) < threshold)
				result[(j*width + i) * 3] = result[(j*width + i) * 3 + 1] = result[(j*width + i) * 3 + 2] = 0;
			else
				result[(j*width + i) * 3] = result[(j*width + i) * 3 + 1] = result[(j*width + i) * 3 + 2] = 255.0;
		}
	}
	delete[] F;
	delete[] gaussC;
	delete[] gaussS;
	delete[] gaussM;
}
static int key[100][2];
static int vis[1000][1000];
static int keyp[1000][1000];
static int num=0,len=0;
void FDoG::DecreaseLine()
{
	for(int i=0;i<len;i++)//由于face++标出的特征点有偏差，在这里做一定程度的8连通扩展，实现面部轮廓清晰的目标
	{
		for(int j=1;j<extend;j++){
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]-j,key[i][1]-j);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]-j,key[i][1]);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]-j,key[i][1]+j);
		memset(vis,0,sizeof(vis));
        kDFS(key[i][0],key[i][1]-j);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0],key[i][1]);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0],key[i][1]+j);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]+j,key[i][1]+j);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]+j,key[i][1]);
		memset(vis,0,sizeof(vis));
		kDFS(key[i][0]+j,key[i][1]-j);
		}
	}
	for(int j = 1;j < height-1;j++)
	{
		for(int i = 1;i < width-1;i++)
			{
				if(result[(j*width + i) * 3]==255||keyp[j][i])continue;//保留面部特征轮廓
				if(result[(j*width + i) * 3]==0)
				{
					memset(vis,0,sizeof(vis));
					DFS(j,i);
					if(num<deleteline)_DFS(j,i);//阈值去线，简化图片
					num=0;
				}
		     }
	
	}
	for(int j = 1;j < height-1;j++)
	{
		for(int i = 1;i < width-1;i++)
			{
				if(result[(j*width + i) * 3]==255)continue;
				if(result[(j*width + i) * 3]==0)
				{
					memset(vis,0,sizeof(vis));
					DFS(j,i);
					if(num<30)_DFS(j,i);//二次阈值去线，简化图片
					num=0;
				}
		     }
	}
}


void FDoG::DFS(int i,int j)
{
	if(num>=deleteline)return;
	if(vis[i][j]==1||i<0||i>=height||j<0||j>=width)return;
	if(result[(i*width + j) * 3]==0)
	{
		 vis[i][j]=1;
		 num++;
		 DFS(i-1,j-1);
		 DFS(i-1,j);
		 DFS(i-1,j+1);
		 DFS(i,j-1);
		 DFS(i,j+1);
		 DFS(i+1,j-1);
		 DFS(i+1,j);
		 DFS(i+1,j+1);
	}  
}

void FDoG::_DFS(int i,int j)//删除线条DFS
{
	if(i<0||i>=height||j<0||j>=width)return;
	if(result[(i*width + j) * 3]==0)
	{
		
		 result[(i*width + j) * 3] = result[(i*width + j) * 3 + 1] = result[(i*width + j) * 3 + 2] = 255.0;
		 _DFS(i-1,j-1);
		 _DFS(i-1,j);
		 _DFS(i-1,j+1);
		 _DFS(i,j-1);
		 _DFS(i,j+1);
		 _DFS(i+1,j-1);
		 _DFS(i+1,j);
		 _DFS(i+1,j+1);
	}  
}

void FDoG::kDFS(int i,int j)//特征点标记DFS
{
	if(vis[i][j]==1||i<0||i>=height||j<0||j>=width)return;
	
	if(result[(i*width + j) * 3]==0)
	{
		 vis[i][j]=1;
		 keyp[i][j]=1;
		 kDFS(i-1,j-1);
		 kDFS(i-1,j);
		 kDFS(i-1,j+1);
		 kDFS(i,j-1);
		 kDFS(i,j+1);
		 kDFS(i+1,j-1);
		 kDFS(i+1,j);
		 kDFS(i+1,j+1);
	}  
}
//计算卷积核的大小，同时也求出在每一点处应乘上的高斯函数值
double* FDoG::MakeGaussians(int* n, double sigma)
{
	int t = 0;
	do
	{
		t++;
	}while(Gaussian(t, 0.0, sigma) > 1e-3);
	*n = t;
	double* gauss = new double[t];
	for(int i = 0; i < t; i++)
		gauss[i] = Gaussian((double)i, 0.0, sigma);
	return gauss;
}

//计算高斯函数值
inline double FDoG::Gaussian(double x, double mean, double sigma)
{
	return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) / sqrt(PI * 2.0 * sigma * sigma) );
}

inline int FDoG::round(double x)
{
	return (int)(x + 0.5);
}

static int n,m;//顶点个数、边数
static int Edge[MAXN][MAXN],lowcost[MAXN],nearvex[MAXN];
static  struct p
{
	int x;
	int y;
	p(){x=0,y=0;}
}xy[MAXN][30000],point[MAXN][MAXN];//xy用于记录各标记连通线条的点集，point用于记录连通线条间边的两个端点

void FDoG::prim(int u0)//最小生成树算法用于形成一笔画
{
    int i,j;
	cout<<"debug   prim"<<endl;
	for(i=1;i<=n;i++)
	{
		lowcost[i]=Edge[u0][i];
		nearvex[i]=u0;
	}
	nearvex[u0]=-1;
	for(i=1;i<n;i++)//树的n-1条边的构造
	{
	   	int min=INF;
		int v=-1;
		for(j=1;j<=n;j++)
		{
			if(nearvex[j]!=-1&&min>lowcost[j])
			{
				v=j; min=lowcost[j];
			}
		}
		if(v!=-1)
		{
          //nearvex[v],v,lowcost[v]
 /////////////////////////bresenham画线/////////////////////////////////////////////////
	  int x0,y0,x1,y1;
	  x0=point[nearvex[v]][v].x,y0=point[nearvex[v]][v].y,x1=point[v][nearvex[v]].x,y1=point[v][nearvex[v]].y;
	  int x,y;
      double d,ke;
	  if(x1-x0==0)
	  {
	   if(y0>y1)
	   {
		  double tx=x0;
          double ty=y0;
		  x0=x1;y0=y1;
		  x1=tx;
		  y1=ty;
	   }
	 x=x0;y=y0;
     d=1;
     for(y=y0;y<=y1;y++)
     {
		 result[(x*width + y) * 3] = result[(x*width + y) * 3 + 1] = result[(x*width + y) * 3 + 2] = 0;
		 for(int j=1;j<10;j++)
		  {
		    result[((x-1)*width + y) * 3] = result[((x-1)*width + y) * 3 + 1] = result[((x-1)*width + y) * 3 + 2] = 0;
			result[(x*width + y-1) * 3] = result[(x*width + y-1) * 3 + 1] = result[(x*width + y-1) * 3 + 2] = 0;
			result[(x*width + y+1) * 3] = result[(x*width + y+1) * 3 + 1] = result[(x*width + y+1) * 3 + 2] = 0;
			result[((x+1)*width + y) * 3] = result[((x+1)*width + y) * 3 + 1] = result[((x+1)*width + y) * 3 + 2] = 0;
		  }
	    }
	  }
	  else 
	  {
		  ke=(y1-y0)/(x1-x0);
	  if(0<=ke && ke<=1)
	 {
       if(x0>x1)
	   {
		  double tx=x0;
          double ty=y0;
		  x0=x1;y0=y1;
		  x1=tx;
		  y1=ty;
	   }
	 x=x0;y=y0;

	 d=0.5-ke;

     for(x=x0;x<=x1;x++)
     {
		 result[(x*width + y) * 3] = result[(x*width + y) * 3 + 1] = result[(x*width + y) * 3 + 2] = 0;
		  for(int j=1;j<10;j++)
		  {
		    result[((x-1)*width + y) * 3] = result[((x-1)*width + y) * 3 + 1] = result[((x-1)*width + y) * 3 + 2] = 0;
			result[(x*width + y-1) * 3] = result[(x*width + y-1) * 3 + 1] = result[(x*width + y-1) * 3 + 2] = 0;
			result[(x*width + y+1) * 3] = result[(x*width + y+1) * 3 + 1] = result[(x*width + y+1) * 3 + 2] = 0;
			result[((x+1)*width + y) * 3] = result[((x+1)*width + y) * 3 + 1] = result[((x+1)*width + y) * 3 + 2] = 0;
		  }
          if(d<0)
		  {
			y++;
          	d+=1-ke;
		  }
		else 
            d-=ke;
      }
     }
     if(ke>1)
	{
       if(y0>y1)
	   {
		  double tx=x0;
          double ty=y0;
		  x0=x1;y0=y1;
		  x1=tx;
		  y1=ty;
	   }
	 x=x0;y=y0;
     d=1-0.5*ke;
     for(y=y0;y<=y1;y++)
     {
		 result[(x*width + y) * 3] = result[(x*width + y) * 3 + 1] = result[(x*width + y) * 3 + 2] = 0;
		 for(int j=1;j<10;j++)
		  {
		    result[((x-1)*width + y) * 3] = result[((x-1)*width + y) * 3 + 1] = result[((x-1)*width + y) * 3 + 2] = 0;
			result[(x*width + y-1) * 3] = result[(x*width + y-1) * 3 + 1] = result[(x*width + y-1) * 3 + 2] = 0;
			result[(x*width + y+1) * 3] = result[(x*width + y+1) * 3 + 1] = result[(x*width + y+1) * 3 + 2] = 0;
			result[((x+1)*width + y) * 3] = result[((x+1)*width + y) * 3 + 1] = result[((x+1)*width + y) * 3 + 2] = 0;
		  }
        if(d>=0)
		{
			x++;
         	d+=1-ke;
		}
		else 
            d+=1;
       }
	 }
     if(ke<-1)
	 {
       if(y0<y1)
	   {
		  double tx=x0;
          double ty=y0;
		  x0=x1;y0=y1;
		  x1=tx;
		  y1=ty;
	   }
	 x=x0;y=y0;
     d=-1-0.5*ke;
     for(y=y0;y>y1;y--)
     {
		result[(x*width + y) * 3] = result[(x*width + y) * 3 + 1] = result[(x*width + y) * 3 + 2] = 0;
		for(int j=1;j<10;j++)
		  {
		    result[((x-1)*width + y) * 3] = result[((x-1)*width + y) * 3 + 1] = result[((x-1)*width + y) * 3 + 2] = 0;
			result[(x*width + y-1) * 3] = result[(x*width + y-1) * 3 + 1] = result[(x*width + y-1) * 3 + 2] = 0;
			result[(x*width + y+1) * 3] = result[(x*width + y+1) * 3 + 1] = result[(x*width + y+1) * 3 + 2] = 0;
			result[((x+1)*width + y) * 3] = result[((x+1)*width + y) * 3 + 1] = result[((x+1)*width + y) * 3 + 2] = 0;
		  }
        if(d<0)
		{
			x++;
         	d-=1+ke;
		}
		else 
            d-=1;
      }
	}
    if(-1<=ke && ke<0)
	{
       if(x0>x1)
	   {
		  double tx=x0;
          double ty=y0;
		  x0=x1;y0=y1;
		  x1=tx;
		  y1=ty;
	   }
	x=x0;y=y0;
    d=-0.5-ke;
    for(x=x0;x<=x1;x++)
    {
		result[(x*width + y) * 3] = result[(x*width + y) * 3 + 1] = result[(x*width + y) * 3 + 2] = 0;
		for(int j=1;j<10;j++)
		  {
		    result[((x-1)*width + y) * 3] = result[((x-1)*width + y) * 3 + 1] = result[((x-1)*width + y) * 3 + 2] = 0;
			result[(x*width + y-1) * 3] = result[(x*width + y-1) * 3 + 1] = result[(x*width + y-1) * 3 + 2] = 0;
			result[(x*width + y+1) * 3] = result[(x*width + y+1) * 3 + 1] = result[(x*width + y+1) * 3 + 2] = 0;
			result[((x+1)*width + y) * 3] = result[((x+1)*width + y) * 3 + 1] = result[((x+1)*width + y) * 3 + 2] = 0;
		  }
        if(d>0)
		{
			y--;
         	d-=1+ke;
		}
		else 
            d-=ke;
     }
   }
}
//////////////////////////////////////////////////////////////////////////////////
			nearvex[v]=-1;
		    for(j=1;j<=n;j++)
		   {
			 if(nearvex[j]!=-1&&lowcost[j]>Edge[v][j])
			 {
				 lowcost[j]=Edge[v][j];
				 nearvex[j]=v;
			 }
		   }
		}
		
	}
}

static int T=1,visit[1000][1000],pos=0,minint=0x7fffffff;
void FDoG::CountEdge()
{
	T=1,pos=0,minint=0x7fffffff;
	memset(visit,0,sizeof(visit));
	for(int j = 1;j < height-1;j++)
	{
		for(int i = 1;i < width-1;i++)
			{
				
				if(result[(j*width + i) * 3]==255.0||visit[j][i]!=0)continue;
				if(result[(j*width + i) * 3]==0&&visit[j][i]==0){DFS_(j,i),T++;cout<<pos<<endl;pos=0;}//连通域标记
				
		    }
	}
	n=T-1;
	memset( Edge, 0, sizeof(Edge) );
	for(int i=1;i<T;i++)//构造图论模型
	{
		//每一个连通域
		cout<<"find edge"<<i<<endl;

		for(int j=i+1;j<T;j++)//剩余连通域	  
		  {
		     for(int k=0;xy[i][k].x;k++)
		    {
	      
			  for(int t=0;xy[j][t].x;t++)
			  {
				  int tmp=(xy[i][k].x-xy[j][t].x)*(xy[i][k].x-xy[j][t].x)+(xy[i][k].y-xy[j][t].y)*(xy[i][k].y-xy[j][t].y);
				  if(minint>tmp)
				  {  minint=tmp;
				     point[i][j].x=xy[i][k].x,point[i][j].y=xy[i][k].y,point[j][i].x=xy[j][t].x,point[j][i].y=xy[j][t].y;
				  }
			  }
			  Edge[i][j]=minint;
			  Edge[j][i]=minint;
		  }
          minint=0x7fffffff;
		}

	 }
	for(int i=1; i<T; i++ ) //对邻接矩阵中其他元素值进行赋值
   {
      for(int j=1; j<T; j++ )
      {
        if( i==j ) Edge[i][j] = 0;
        else if( Edge[i][j]==0 ) Edge[i][j] = INF;
      }
   }
}
void FDoG::DFS_(int i,int j)
{
	if(visit[i][j]!=0||i<0||i>=height||j<0||j>=width)return;
	if(result[(i*width + j) * 3]==0) 
	{
		visit[i][j]=T;
		 xy[T][pos].x=i;
		 xy[T][pos].y=j;
		 pos++;
		 DFS_(i-1,j-1);
		 DFS_(i-1,j);
		 DFS_(i-1,j+1);
		 DFS_(i,j-1);
		 DFS_(i,j+1);
		 DFS_(i+1,j-1);
		 DFS_(i+1,j);
		 DFS_(i+1,j+1);
	}
}