//#include <iostream>
//#include <opencv2/opencv.hpp>
//#include <opencv2/core/types_c.h>
//
//using namespace cv;
//using namespace std;
//
//Mat fish, fishROI, nor; //鱼眼图（原图），ROI区域，正常图（效果图）
//Mat mx, my;            //remap参数，变换矩阵
//Mat Rotation;         //旋转后的坐标矩阵
//Mat rou, phi;          //rou 和 phi
//
//int ox, oy;            //原点xy坐标
//const int num_View = 6;
//
//double f = 1, DR = 180, dr = 20, thetaMax = 80, thetaMin = 10, r = 180;
//double viewPoint = 45;                               //虚拟视点度数
//double view_phi[num_View + 1];                       //分割区域的度数
//double points[num_View + 1];                         //虚拟视点的投影方向
//
//void setO(Mat& m)
//{
//	ox = m.cols / 2;
//	oy = m.rows / 2;
//}
//
//
//
//int main()
//{
//	fish = imread("1.jpg", 1);
//	nor.create(fish.size(), fish.type());
//	mx.create(fish.size(), CV_32FC1);
//	my.create(fish.size(), CV_32FC1);
//	setO(fish);
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	vector<Point2f> point;
//	//i对应x  j对应y
//	for (int j = 0; j<fish.rows; j++)
//	{
//		for (int i = 0; i<fish.cols; i++)
//		{
//			point.push_back(Point2f(i + 1 - ox, oy - j + 1));
//			//cout<<i+1-ox<<","<<oy-j+1<<"     "<<i<<","<<j<<endl;
//		}
//	}
//
//	//极坐标与直角坐标互相转换的xy数组
//	Mat xp(point.size(), 1, CV_32F, &point[0].x, 2 * sizeof(float));
//	Mat yp(point.size(), 1, CV_32F, &point[0].y, 2 * sizeof(float));
//
//	//极坐标rou与phi
//	cartToPolar(xp, yp, rou, phi);
//	//cout<<phi<<endl;
//	//ps:角度转弧度 π/180×角度,弧度变角度 180/π×弧度
//	//ps:依行遍历
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	f = 147;
//	for (int j = 0; j < nor.rows; j++)
//	{
//		for (int i = 0; i < nor.cols; i++)
//		{
//			double range = rou.at<float>(i + j*fish.cols);
//			mx.at<float>(j, i) = static_cast<float>((ox + ((i - ox) / ((range / f) / atan(range / f)))));
//			my.at<float>(j, i) = static_cast<float>((oy + ((j - oy) / ((range / f) / atan(range / f)))));
//		}
//	}
//	remap(fish, nor, mx, my, CV_INTER_LINEAR);
//	imshow("1", nor);
//	imshow("2", fish);
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	//每个区域的度数
//	for (int i = 0; i<num_View; i++)
//	{
//		view_phi[i] = i * 2 * CV_PI / num_View;//弧度制
//		//cout<<view_phi[i]*180/CV_PI<<endl;
//	}
//	view_phi[num_View] = view_phi[0] + 2 * CV_PI;
//	//cout<<view_phi[num_View]*180/CV_PI<<endl;
//
//	//视点的位置角度
//	for (int i = 0; i<num_View; i++)
//	{
//		points[i] = ((view_phi[i + 1] + view_phi[i]) / 2) * 180 / CV_PI;//角度制
//		//cout<<points[i]<<endl;
//	}
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	//极坐标映射到球坐标：二维->三维   ***
//	Mat xx(point.size(), 1, CV_32F);
//	Mat yy(point.size(), 1, CV_32F);
//	Mat zz(point.size(), 1, CV_32F);
//	double *tphi = new double[point.size()];
//	double *ttheta = new double[point.size()];
//	for (int i = 0; i<point.size(); i++)
//	{
//		//        ttheta[i]=(CV_PI/2-rou.at<float>(i)/f);
//		//        tphi[i]=(CV_PI+phi.at<float>(i));
//
//		//处理ROI
//		if (rou.at<float>(i)<190)
//		{
//			ttheta[i] = (CV_PI / 2 - rou.at<float>(i) / f);
//			tphi[i] = (phi.at<float>(i)-CV_PI);
//
//		}
//		else
//		{
//			tphi[i] = 0;
//			ttheta[i] = 0;
//		}
//		//cout<<rou.at<float>(i)<<endl;
//	}
//	for (int j = 0; j<fish.rows; j++)
//	{
//		for (int i = 0; i<fish.cols; i++)
//		{
//			xx.at<float>(i + j*fish.cols) = static_cast<float>(r*cos(ttheta[i + j*fish.cols])*-cos(tphi[i + j*fish.cols]));
//			yy.at<float>(i + j*fish.cols) = static_cast<float>(r*cos(ttheta[i + j*fish.cols])*sin(tphi[i + j*fish.cols]));
//			zz.at<float>(i + j*fish.cols) = static_cast<float>(r*sin(ttheta[i + j*fish.cols]));
//			//cout<<xx.at<float>(i+j*fish.cols)<<"a= "<<yy.at<float>(i+j*fish.cols)<<"b= "<<zz.at<float>(i+j*fish.cols)<<"     "<<i<<","<<j<<endl;
//		}
//	}
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	//[x,y,z]
//	Mat xyz(point.size(), 3, CV_32F);
//	for (int i = 0; i<point.size(); i++)
//	{
//		xyz.at<float>(i * 3) = xx.at<float>(i);
//		xyz.at<float>(i * 3 + 1) = yy.at<float>(i);
//		xyz.at<float>(i * 3 + 2) = zz.at<float>(i);
//	}
//	
//	//旋转矩阵
//	double mtheta = CV_PI / 4, mphi = -CV_PI / 4;
//	double r1[3][3] = { cos(mtheta), sin(mtheta), 0, -sin(mtheta), cos(mtheta), 0, 0, 0, 1 };
//	double r2[3][3] = { cos(mphi), 0, -sin(mphi), 0, 1, 0, sin(mphi), 0, cos(mphi) };
//	double r3[3][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, };
//	double **db_xyz = new double *[point.size()];
//	double **db_xyz_2 = new double *[point.size()];
//	for (int i = 0; i<point.size(); i++)
//	{
//		db_xyz[i] = new double[3];
//		db_xyz_2[i] = new double[3];
//		db_xyz[i][0] = xyz.at<float>(i * 3);
//		db_xyz[i][1] = xyz.at<float>(i * 3 + 1);
//		db_xyz[i][2] = xyz.at<float>(i * 3 + 2);
//		//cout<<db_xyz[i][0]<<"a= "<<db_xyz[i][1]<<"b= "<<db_xyz[i][2]<<endl;
//	}
//	//    CvMat R1=cvMat(3,3,CV_32F,r1);
//	//    CvMat R2=cvMat(3,3,CV_32F,r2);
//	//    CvMat R3=cvMat(3,3,CV_32F,r3);
//	//    CvMat XYZ=cvMat(point.size(),3,CV_32F,db_xyz);
//	//    CvMat XYZ_2=cvMat(point.size(),3,CV_32F,db_xyz_2);
//	for (int i = 0; i<3; i++)
//	{
//		r3[i][0] = r1[i][0] * r2[0][0] + r1[i][1] * r2[1][0] + r1[i][2] * r2[2][0];
//		r3[i][1] = r1[i][0] * r2[0][1] + r1[i][1] * r2[1][1] + r1[i][2] * r2[2][1];
//		r3[i][2] = r1[i][0] * r2[0][2] + r1[i][1] * r2[1][2] + r1[i][2] * r2[2][2];
//	}
//	//cout<<"b"<<endl;
//	for (int i = 0; i<point.size(); i++)
//	{
//		db_xyz_2[i][0] = db_xyz[i][0] * r3[0][0] + db_xyz[i][1] * r3[1][0] + db_xyz[i][2] * r3[2][0];
//		db_xyz_2[i][1] = db_xyz[i][0] * r3[0][1] + db_xyz[i][1] * r3[1][1] + db_xyz[i][2] * r3[2][1];
//		db_xyz_2[i][2] = db_xyz[i][0] * r3[0][2] + db_xyz[i][1] * r3[1][2] + db_xyz[i][2] * r3[2][2];
//	}
//	//cout<<"b"<<endl;
//	for (int i = 0; i<point.size(); i++)
//	{
//		xx.at<float>(i) = static_cast<float>(db_xyz_2[i][0]);
//		yy.at<float>(i) = static_cast<float>(db_xyz_2[i][1]);
//		zz.at<float>(i) = static_cast<float>(db_xyz_2[i][2]);
//		//cout<<xx.at<float>(i)<<"a= "<<yy.at<float>(i)<<"b= "<<zz.at<float>(i)<<endl;
//	}
//	/*--------------------------------------------------------------------------------------------------------------------*/
//	//球坐标映射到平面坐标：三维->二维
//	Mat nrou(point.size(), 1, CV_32F);
//	Mat nphi(point.size(), 1, CV_32F);
//	for (int j = 0; j<fish.rows; j++)
//	{
//		for (int i = 0; i<fish.cols; i++)
//		{
//			nrou.at<float>(i + j*fish.cols) = static_cast<float>(((CV_PI / 2 - atan((zz.at<float>(i + j*fish.cols) / (sqrt(xx.at<float>(i + fish.cols)*xx.at<float>(i + fish.cols) + yy.at<float>(i + fish.cols)*yy.at<float>(i + fish.cols))))))*f));
//			nphi.at<float>(i + j*fish.cols) = static_cast<float>(CV_PI + atan2(yy.at<float>(i + j*fish.cols), xx.at<float>(i + j*fish.cols)));
//			//cout<<xx.at<float>(i+j*fish.cols)<<"a= "<<yy.at<float>(i+j*fish.cols)<<"b= "<<zz.at<float>(i+j*fish.cols)<<" "<<i<<" "<<j<<endl;
//			//cout<<(float)((CV_PI/2-atan((zz.at<float>(i+j*fish.cols)/(sqrt(xx.at<float>(i+fish.cols)*xx.at<float>(i+fish.cols)+yy.at<float>(i+fish.cols)*yy.at<float>(i+fish.cols))))))*f)<<endl;
//			//cout<<CV_PI+atan(yy.at<float>(i+j*fish.cols)/xx.at<float>(i+j*fish.cols))<<" "<<endl;
//			//cout<<nrou.at<float>(i+j*fish.cols)<<"  "<<nphi.at<float>(i+j*fish.cols)<<endl;
//
//			mx.at<float>(j, i) = static_cast<float>(xx.at<float>(i + j*fish.cols) + ox);
//			my.at<float>(j, i) = static_cast<float>(yy.at<float>(i + j*fish.cols) + oy);
//		}
//	}
//	remap(fish, nor, mx, my, CV_INTER_LINEAR);
//	imshow("3", nor);
//
//	for (int j = 0; j<fish.rows; j++)
//	{
//		for (int i = 0; i<fish.cols; i++)
//		{
//			double range = rou.at<float>(i + j*fish.cols);
//			//            mx.at<float>(j,i) = static_cast<float>(ox+(i-ox)/((nrou.at<float>(i+j*fish.cols)/f)));
//			//            my.at<float>(j,i) = static_cast<float>(oy+(j-oy)/((nrou.at<float>(i+j*fish.cols)/f)));
//			mx.at<float>(j, i) = static_cast<float>((ox + ((i - ox) / ((range / f) / atan(range / f)))));
//			my.at<float>(j, i) = static_cast<float>((oy + ((j - oy) / ((range / f) / atan(range / f)))));
//			//            cout<<phi.at<float>(i+j*fish.cols)<<"  "<<nphi.at<float>(i+j*fish.cols)<<endl;
//			//            cout<<((range/f)/atan(range/f))<<"  "<<tan(nrou.at<float>(i+j*fish.cols)/f)<<endl;
//		}
//	}
//	remap(nor, fishROI, mx, my, CV_INTER_LINEAR);
//	imshow("4", fishROI);
//
//	delete[](tphi);
//	delete[](ttheta);
//	delete[](db_xyz);
//	delete[](db_xyz_2);
//	cvWaitKey();
//	return 0;
//}