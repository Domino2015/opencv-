//#include<opencv2\opencv.hpp>
//#include<iostream>
//#include<math.h>
//
//using namespace std;
//using namespace cv;
//
//Mat srcImg;//ԭͼƬ
//Mat srcImagROI;//ѡȡ����Ȥ����
//Mat outImg;//����Ч��ͼ
//Mat map_x, map_y;
//int ox, oy;//ͼ������
//double phi; //������Ƕ�
//double range;//������뾶
//double rangeMax;//���������뾶
//double newRange;//��һ����r
//int f;//����
//int dr=45;//���ĺ�Բ�뾶
//int DR=200;//��߽��Բ�뾶
//double theta;//�Ƕ�
//double thetamax = 80 * CV_PI / 180;//���Ƕ�
//double thetamin;//��С�Ƕ�
//double viewpoint;//��������ӵ�
//
//int main()
//{
//	//����ԭʼͼ
//	Mat srcImg2 = imread("1111.jpg");
//	srcImg = imread("1.jpg");
//	if (!srcImg.data) { cout << "��ȡͼƬ������ȷ��Ŀ¼���Ƿ���imread����ָ����ͼƬ����~�� \n"; return -1; }
//	imshow("ԭʼͼ", srcImg);
//	//����ͼƬ����Ϊԭ��
//	oy = srcImg.cols / 2;
//	ox = srcImg.rows / 2;
//	//������ԭʼͼһ����Ч��ͼ��x��ӳ��ͼ��y��ӳ��ͼ
//	srcImg.create(srcImg.size(), srcImg.type());
//	map_x.create(srcImg.size() , CV_32FC1);
//	map_y.create(srcImg.size() , CV_32FC1);
//
//	//˫��ѭ��������ÿһ�����ص㣬�ı�map_x & map_y��ֵ
//	for (int j = 1; j < srcImg.rows; j++)
//	{
//		for (int i = 1; i < srcImg.cols; i++)
//		{
//			if ((i - ox) != 0 && (i - oy) != 0)
//			{
//			//����ϵ ԭ��ƽ��
//			double newX = i - ox;
//			double newY = j - oy;
//
//			//ֱ������ת��Ϊ������
//			phi = atan2(newY, newX);//������Ƕ�
//			range = hypot(newX, newY);
//		
//			////�󽹾�
//			//f = DR / thetamax;
//			////�Ƕ�
//			//theta = range / f;
//			//thetamin = dr / f;
//			////ģ��һ�������ӽ�
//			//viewpoint = 0.5*(thetamax + thetamin); //��������ӵ�
//			double points = CV_PI;//�����ӽ�
//			double  theta = points/4;//0��
//
//			//������ת��Ϊ��������
//			double qx = range*cos(phi);
//			double qy = range*sin(phi);
//			double qz = range*cos(theta);
//			//�仯
//			double nx = qz*sin(theta) + qx*cos(theta);
//
//
//			//double range_1 = sqrt(qx*qx+qz*qz);//�任����µļ�����İ뾶
//			//double theta_1 = CV_PI / 2 - atan(qy / range);
//			// double phi_2 = acos(qx / range);//�任����µļ�����ĽǶ�
//			map_x.at<float>(j, i) = static_cast<float>(nx+ox);
//			map_y.at<float>(j, i) = static_cast<float>(qy + oy);
//			//�ı�map_x & map_y��ֵ. 
//			// map_x.at<float>(j, i) = static_cast<float>(range*cos(phi) + ox);
//			 //map_y.at<float>(j, i) = static_cast<float>(range*sin(phi) + oy);
//
//
//		}
//	}
//	//������ӳ�����
//		remap(srcImg, outImg, map_x, map_y, CV_INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0));
//
//	//��ʾЧ��ͼ
//	imshow("�����򴰿ڡ�", outImg);
//	} //ifend
//	waitKey();
//}