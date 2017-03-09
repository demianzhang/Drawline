环境VS2010，OpenCV2.1.0
需要matlab调用face++的接口得到的key.txt文件

main函数里的参数初始化：
针对不同图片要选择合适的最后两个参数进行初始化

注：这里一直没有找到bug，
对fdog.CountEdge();fdog.prim(1);
进行多次迭代也没有完全连通，效果不是特别好
理论上会完全连通的，Audrey Hepburn.jpg的连通效果比lena.jpg要好很多
