//#include<opencv2\opencv.hpp>
//#include<iostream>
//#include<math.h>
//
//using namespace std;
//using namespace cv;
//
//Mat srcImg;//原图片
//Mat srcImagROI;//选取感兴趣区域
//Mat outImg;//最终效果图
//Mat map_x, map_y;
//int ox, oy;//图像中心
//double phi; //极坐标角度
//double range;//极坐标半径
//double rangeMax;//极坐标最大半径
//double newRange;//归一化后r
//int f;//焦距
//int dr=45;//中心黑圆半径
//int DR=200;//外边界大圆半径
//double theta;//角度
//double thetamax = 80 * CV_PI / 180;//最大角度
//double thetamin;//最小角度
//double viewpoint;//多个虚拟视点
//
//int main()
//{
//	//载入原始图
//	Mat srcImg2 = imread("1111.jpg");
//	srcImg = imread("1.jpg");
//	if (!srcImg.data) { cout << "读取图片错误，请确定目录下是否有imread函数指定的图片存在~！ \n"; return -1; }
//	imshow("原始图", srcImg);
//	//设置图片中心为原点
//	oy = srcImg.cols / 2;
//	ox = srcImg.rows / 2;
//	//创建和原始图一样的效果图，x重映射图，y重映射图
//	srcImg.create(srcImg.size(), srcImg.type());
//	map_x.create(srcImg.size() , CV_32FC1);
//	map_y.create(srcImg.size() , CV_32FC1);
//
//	//双层循环，遍历每一个像素点，改变map_x & map_y的值
//	for (int j = 1; j < srcImg.rows; j++)
//	{
//		for (int i = 1; i < srcImg.cols; i++)
//		{
//			if ((i - ox) != 0 && (i - oy) != 0)
//			{
//			//坐标系 原点平移
//			double newX = i - ox;
//			double newY = j - oy;
//
//			//直角坐标转换为极坐标
//			phi = atan2(newY, newX);//极坐标角度
//			range = hypot(newX, newY);
//		
//			////求焦距
//			//f = DR / thetamax;
//			////角度
//			//theta = range / f;
//			//thetamin = dr / f;
//			////模拟一个虚拟视角
//			//viewpoint = 0.5*(thetamax + thetamin); //多个虚拟视点
//			double points = CV_PI;//鱼眼视角
//			double  theta = points/4;//0角
//
//			//极坐标转换为球面坐标
//			double qx = range*cos(phi);
//			double qy = range*sin(phi);
//			double qz = range*cos(theta);
//			//变化
//			double nx = qz*sin(theta) + qx*cos(theta);
//
//
//			//double range_1 = sqrt(qx*qx+qz*qz);//变换后的新的极坐标的半径
//			//double theta_1 = CV_PI / 2 - atan(qy / range);
//			// double phi_2 = acos(qx / range);//变换后的新的极坐标的角度
//			map_x.at<float>(j, i) = static_cast<float>(nx+ox);
//			map_y.at<float>(j, i) = static_cast<float>(qy + oy);
//			//改变map_x & map_y的值. 
//			// map_x.at<float>(j, i) = static_cast<float>(range*cos(phi) + ox);
//			 //map_y.at<float>(j, i) = static_cast<float>(range*sin(phi) + oy);
//
//
//		}
//	}
//	//进行重映射操作
//		remap(srcImg, outImg, map_x, map_y, CV_INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0));
//
//	//显示效果图
//	imshow("【程序窗口】", outImg);
//	} //ifend
//	waitKey();
//}