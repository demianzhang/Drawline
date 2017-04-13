# Portrait Printing Algorithm and Optimization 
肖像打印算法设计与优化

功能：人物肖像轮廓提取，规划线条路径，动态展示

## Intro
| File | Description | Dependencies |
| --- | --- | --- | 
| **ETF+FDoG** | Coherent Line Drawing 算法提取照片轮廓，简化线条，连通域互连 | `OpenCV2.1.0` |
| **Face++_Matlab** | 脸部关键点提取 | `matlab` |
| **Drawline** | 线条细化抽取骨架，DFS动态打印 | `OpenCV2.4.9` |

## Run
```bash
Face++_Matlab ---> ETF+FDoG ---> Drawline
```
需要matlab调用face++的接口得到的key.txt文件  
ETF+FDoG的main函数里的参数初始化：针对不同图片要选择合适的最后两个参数进行初始化  
选择不同的图片进行绘制时，图片源为ETF+FDoG程序的运行结果(二值图)  
同时要自己添加一个与原图片大小相同的白色背景图片background.jpg


Face++ Research Toolkit - Matlab SDK      
Copyright © 2013 Megvii, Inc. All Rights Reserved.   
