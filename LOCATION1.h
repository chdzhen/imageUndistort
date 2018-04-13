
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <time.h>
#include <windows.h>

#include "test.h"


//double to string
string doubleToString(double num)
{
	/*char str[256];
	sprintf(str, "%.0f", num);
	string result = str;
	return result;*/
	return to_string(int(num));

}


// location the cooder
void getcoorder(cv::Mat oriRed, vector<cv::Point>& colcooder)
{
	vector<vector<cv::Point>> contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(oriRed, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
	for (int i = 0; i < contours.size(); i++)
	{
		if (contours[i].size() < 5)
		{
			continue;
		}
		cv::RotatedRect box = cv::fitEllipse(contours[i]);
		colcooder.push_back(box.center);
	}
}


//get center by num of contours
cv::Point getcenter(cv::Mat image)
{
	vector<vector<cv::Point>> contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(image, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0));
	int temp = 0;
	int maxindex = 0;
	for (int i = 0; i < contours.size(); i++)
	{
		if (contours[i].size() > temp)
		{
			temp = contours[i].size();
			maxindex = i;
		}
	}
	cv::Point center;

	if (contours[maxindex].size() >= 5)
	{
		cv::RotatedRect box = cv::fitEllipse(contours[maxindex]);
		center = box.center;
	}
	return center;
}


// sort potions
void colsort(vector<cv::Point> &colcooder, int mood)
{
	for (int i = 0; i < colcooder.size(); i++)
	{
		for (int j = i; j < colcooder.size(); j++)
		{
			if (mood == 1)
			{
				if (colcooder[i].y > colcooder[j].y)
				{
					cv::Point temp = colcooder[i];
					colcooder[i] = colcooder[j];
					colcooder[j] = temp;
				}
			}
			if (mood == 0)
			{
				if (colcooder[i].x > colcooder[j].x)
				{
					cv::Point temp = colcooder[i];
					colcooder[i] = colcooder[j];
					colcooder[j] = temp;
				}
			}
		}
	}
}


//get margin of positions(row and col)
vector<double> getMargin(vector<cv::Point> colcooder, cv::Point center, int mode)
{
	vector<double> rowMargin;

	int centerNum=0;
	for (int i = 0; i < colcooder.size(); i++)
	{
		if (mode == 0)
		{
			if (abs(colcooder[i].x - center.x)<10)
			{
				centerNum = i;
			}
		}
		if (mode == 1)
		{
			if (abs(colcooder[i].y - center.y)<10)
			{
				centerNum = i;
			}
		}

	}
	for (int ii = 1; ii < centerNum + 1; ii++)
	{
		double temp;
		if (mode == 0)
		{
			temp = colcooder[ii].x - colcooder[ii - 1].x;
		}
		if (mode == 1)
		{
			temp = temp = colcooder[ii].y - colcooder[ii - 1].y;
		}

		rowMargin.push_back(temp);
	}
	return rowMargin;
}




// git fitdata(fitimage)
vector<cv::Point> getFitData(vector<cv::Point> allCooder, cv::Point center, cv::Point fitcenter)
{
	vector<cv::Point> desFit;
	for (int i = 0; i < allCooder.size(); i++)
	{
		double fenzi_x = (allCooder[i].x - center.x);
		double fenmu_x = (center.x + 100);
		double x = fenzi_x / fenmu_x;

		double fenzi_y = (allCooder[i].y - center.y);
		double fenmu_y = (center.y + 100);
		double y = fenzi_y / fenmu_y;

		cv::Point temp;
		temp.x = fitcenter.x + xx::get_answer(x, y)*(center.y + 100);
		temp.y = fitcenter.y + yy::get_answer(x, y)*(center.y + 100);
		//cout << xx::get_answer() << "," << yy::get_answer() << endl;

		desFit.push_back(temp);
	}
	return desFit;
}


void writefitImage(cv::Mat image, vector<cv::Point> fitdata, int mode)
{
	for (int i = 0; i < fitdata.size(); i++)
	{

		if (mode == 0)
		{
			cv::circle(image, fitdata[i], 9, cv::Scalar(0, 0, 255), -1, 8);
		}
		if (mode == 1)
		{
			cv::circle(image, fitdata[i], 9, cv::Scalar(255, 0, 0), -1, 8);
		}
		if (mode == 2)
		{
			cv::circle(image, fitdata[i], 9, cv::Scalar(0, 255, 0), -1, 8);
		}
	}

}

// undistortPoints
void UndistortPoints(cv::InputArray distorted, cv::OutputArray undistorted,
	cv::InputArray K, cv::InputArray D, cv::InputArray R=cv::noArray())
{
	// will support only 2-channel data now for points
	CV_Assert(distorted.type() == CV_32FC2 || distorted.type() == CV_64FC2);
	undistorted.create(distorted.size(), distorted.type());

	//CV_Assert(P.empty() || P.size() == cv::Size(3, 3) || P.size() == cv::Size(4, 3));
	CV_Assert(R.empty() || R.size() == cv::Size(3, 3) || R.total() * R.channels() == 3);
	CV_Assert(D.total() == 4 && K.size() == cv::Size(3, 3) && (K.depth() == CV_32F || K.depth() == CV_64F));

	cv::Vec2d f, c;
	if (K.depth() == CV_32F)
	{
		cv::Matx33f camMat = K.getMat();
		f = cv::Vec2f(camMat(0, 0), camMat(1, 1));
		c = cv::Vec2f(camMat(0, 2), camMat(1, 2));
	}
	else
	{
		cv::Matx33d camMat = K.getMat();
		f = cv::Vec2d(camMat(0, 0), camMat(1, 1));
		c = cv::Vec2d(camMat(0, 2), camMat(1, 2));
	}

	cv::Vec4d k = D.depth() == CV_32F ? (cv::Vec4d)*D.getMat().ptr<cv::Vec4f>() : *D.getMat().ptr<cv::Vec4d>();

	cv::Matx33d RR = cv::Matx33d::eye();
	if (!R.empty() && R.total() * R.channels() == 3)
	{
		cv::Vec3d rvec;
		R.getMat().convertTo(rvec, CV_64F);
		RR = cv::Affine3d(rvec).rotation();
	}
	else if (!R.empty() && R.size() == cv::Size(3, 3))
		R.getMat().convertTo(RR, CV_64F);

	// start undistorting
	const cv::Vec2f* srcf = distorted.getMat().ptr<cv::Vec2f>();
	const cv::Vec2d* srcd = distorted.getMat().ptr<cv::Vec2d>();
	cv::Vec2f* dstf = undistorted.getMat().ptr<cv::Vec2f>();
	cv::Vec2d* dstd = undistorted.getMat().ptr<cv::Vec2d>();

	size_t n = distorted.total();
	int sdepth = distorted.depth();

	for (size_t i = 0; i < n; i++)
	{
		cv::Vec2d pi = sdepth == CV_32F ? (cv::Vec2d)srcf[i] : srcd[i];  // image point
		cv::Vec2d pw((pi[0] - c[0]) / f[0], (pi[1] - c[1]) / f[1]);      // world point

		double scale = 1.0;

		double theta_d = sqrt(pw[0] * pw[0] + pw[1] * pw[1]);
		if (theta_d > 1e-8)
		{
			// compensate distortion iteratively
			double theta = theta_d;
			for (int j = 0; j < 10; j++)
			{
				double theta2 = theta*theta, theta4 = theta2*theta2, theta6 = theta4*theta2, theta8 = theta6*theta2;
				theta = theta_d / (1 + k[0] * theta2 + k[1] * theta4 + k[2] * theta6 + k[3] * theta8);
			}

			scale = std::tan(theta) / theta_d;
		}

		cv::Vec2d pu = pw * scale; //undistorted point

		// reproject
		cv::Vec3d pr = RR * cv::Vec3d(pu[0], pu[1], 1.0); // rotated point optionally multiplied by new camera matrix
		cv::Vec2d fi(pr[0] / pr[2], pr[1] / pr[2]);       // final

		if (sdepth == CV_32F)
			dstf[i] = fi;
		else
			dstd[i] = fi;
	}
}

void myUndistortPoints(std::vector<cv::Point> distorted, std::vector<cv::Point>& undistorted, cv::Mat K, cv::Mat D)
{
	std::vector<cv::Point2f> tmpDistorted, tmpUndistorted;
	for (auto p : distorted)
		//for (int i = 0; i < distorted.size(); i++)
	{
		tmpDistorted.push_back(cv::Point2f((float)p.x, (float)p.y));
	}
	if (tmpDistorted.size() != distorted.size())
		std::cout << "point size error!" << std::endl;

	UndistortPoints(tmpDistorted, tmpUndistorted, K, D);

	for (auto p : tmpUndistorted)
	{
		//cout << p << endl;
		//double tt = p.x * K.at<double>(0, 0) + K.at<double>(0, 2);
		//double tr = p.y*K.at<double>(1, 1) + K.at<double>(1, 2);
		undistorted.push_back(cv::Point(p.x * K.at<double>(0, 0) + K.at<double>(0, 2),
			p.y*K.at<double>(1, 1) + K.at<double>(1, 2)));
	}
	if (undistorted.size() != tmpUndistorted.size())
		std::cout << "out point size Error." << std::endl;
}



// add the contrast of image
void ContrastImage(cv::Mat &image, int alpha, int beta)
{
	for (int c = 0; c < 3; c++)
	{
		for (int i = 0; i < image.rows; i++)
		{
			for (int j = 0; j < image.cols; j++)
			{
				int temp;
				if (c == 1)
				{
					temp = (image.at<cv::Vec3b>(i, j)[c] + beta - 50)*alpha;
				}
				else
				{
					temp = (image.at<cv::Vec3b>(i, j)[c] + beta)*alpha;
				}

				if (temp>255)
				{
					temp = 255;
				}
				else if (temp < 0)
				{
					temp = 0;
				}
				image.at<cv::Vec3b>(i, j)[c] = temp;
			}
		}
	}

}

//cameracooder RIGHT
vector<cv::Point> cameraCooderRight(vector<cv::Point>camera_Redcooder, cv::Point cameraCenter, int mode)
{
	vector<cv::Point> camera_Redcooder_right;
	for (int i = 0; i < camera_Redcooder.size(); i++)
	{
		if (mode == 1)
		{
			if (abs(camera_Redcooder[i].x - cameraCenter.x) < 50)
			{
				camera_Redcooder_right.push_back(camera_Redcooder[i]);
			}
		}
		if (mode == 0)
		{
			if (abs(camera_Redcooder[i].y - cameraCenter.y) < 50)
			{
				camera_Redcooder_right.push_back(camera_Redcooder[i]);
			}
		}

	}
	return camera_Redcooder_right;
}


vector<cv::Point> otherCooderRight(const vector<cv::Point> camera_Redcooder_right,
	const vector<cv::Point> camera_Bluecooder_right, vector<cv::Point> camera_Greencooder)
{
	cout << camera_Greencooder.size() << endl;

	for (int i = 0; i < camera_Greencooder.size(); i++)
	{
		for (int j = 0; j < camera_Redcooder_right.size(); j++)
		{
			double lenght1 = sqrt(pow(camera_Greencooder[i].x - camera_Redcooder_right[j].x, 2) +
				pow(camera_Greencooder[i].y - camera_Redcooder_right[j].y, 2));
			if (lenght1 < 50)
			{
				camera_Greencooder.erase(camera_Greencooder.begin() + i);
				i--;
			}

		}

	}
	cout << camera_Greencooder.size() << endl;
	for (int ii = 0; ii < camera_Greencooder.size(); ii++)
	{
		for (int k = 0; k < camera_Bluecooder_right.size(); k++)
		{
			double lenght2 = sqrt(pow(camera_Greencooder[ii].x - camera_Bluecooder_right[k].x, 2) +
				pow(camera_Greencooder[ii].y - camera_Bluecooder_right[k].y, 2));
			if (lenght2 < 50)
			{
				camera_Greencooder.erase(camera_Greencooder.begin() + ii);
				ii--;
			}
		}
	}
	cout << camera_Greencooder.size() << endl;
	return camera_Greencooder;
}


cv::Point ChangePoint(cv::Point beginPoint, int x_length, int y_length)
{
	cv::Point afterPoint;
	afterPoint = cv::Point(beginPoint.x - x_length, beginPoint.y - y_length);
	return afterPoint;
}


double getAverage(vector<double> Margin)
{
	double sum = 0;
	for (int i = 0; i < Margin.size(); i++)
	{
		sum += Margin[i];
	}

	double average = 0;
	average = sum / Margin.size();
	return average;
}


vector<vector<double>> normalizeImageCooder(vector<cv::Point> cameraother_JZ, cv::Point cameracenter_JZ, double rowaverageJZ, double colaverageJZ)
{
	vector<vector<double>> normalizeJZ;


	for (int i = 0; i < cameraother_JZ.size(); i++)
	{
		double temp_x = (cameraother_JZ[i].x - cameracenter_JZ.x) / rowaverageJZ;
		double temp_y = (cameraother_JZ[i].y - cameracenter_JZ.y) / colaverageJZ;

		int rowJZtemp;
		int colJZtemp;

		if (temp_x >= 0)
		{
			rowJZtemp = int(temp_x + 0.5);
		}
		else
		{
			rowJZtemp = int(temp_x - 0.5);
		}

		if (temp_y >= 0)
		{
			colJZtemp = int(temp_y + 0.5);
		}
		else
		{
			colJZtemp = int(temp_y - 0.5);
		}
		vector<double> temp;
		temp.push_back(cameraother_JZ[i].x);
		temp.push_back(cameraother_JZ[i].y);
		temp.push_back(rowJZtemp);
		temp.push_back(colJZtemp);

		normalizeJZ.push_back(temp);
	}

	return normalizeJZ;
}

// calclate the cooder of rotate
vector<cv::Point> CalrotatePoint(vector<cv::Point> points, cv::Point center, double angle)
{
	vector<cv::Point> revector;

	for (int i = 0; i < points.size(); i++)
	{
		double R = sqrt(pow(points[i].x - center.x, 2) + pow(points[i].y - center.y, 2));
		cv::Point tempPoint;
		double rotateangle = atan2(points[i].y - center.y, points[i].x - center.x) - angle;
		tempPoint.y = R*sin(rotateangle) + center.y;
		tempPoint.x = R*cos(rotateangle) + center.x;
		revector.push_back(tempPoint);
	}
	return revector;
}
//
void  writelineandtext(string filename, cv::Mat fitImage, vector<vector<double>> normalizeFit, vector<cv::Point> desCornersFit_row, vector<cv::Point> desCornersFit_col)
{
	for (int i = 0; i < desCornersFit_row.size(); i++)
	{
		cv::line(fitImage, cv::Point(desCornersFit_row[i].x, 0), cv::Point(desCornersFit_row[i].x, fitImage.rows), cv::Scalar(255, 255, 255), 1, 8);
	}

	for (int i = 0; i < desCornersFit_col.size(); i++)
	{
		cv::line(fitImage, cv::Point(0, desCornersFit_col[i].y), cv::Point(fitImage.cols, desCornersFit_col[i].y), cv::Scalar(255, 255, 255), 1, 8);
	}

	int font_face = cv::FONT_HERSHEY_COMPLEX;
	for (int i = 0; i < normalizeFit.size(); i++)
	{
		string text = doubleToString(normalizeFit[i][2]).append(",").append(doubleToString(normalizeFit[i][3]));
		cv::putText(fitImage, text, cv::Point(normalizeFit[i][0] - 5, normalizeFit[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	//cv::imshow("test", fitImage);
	cv::imwrite(filename.append(".jpg"), fitImage);
}



/*************************************** test *************************************************/
//分离四部分
vector<vector<cv::Point>> separatefour(vector<cv::Point> other, cv::Point center)
{
	vector<vector<cv::Point>> refourvector;
	vector<cv::Point> left_up, left_down, right_up, right_down;
	for (int i = 0; i<other.size(); i++)
	{
		//left_up
		if (other[i].x<center.x && other[i].y<center.y)
		{
			left_up.push_back(other[i]);
		}
		//left_dowm
		if (other[i].x>center.x && other[i].y<center.y)
		{
			right_up.push_back(other[i]);
		}
		//right_up
		if (other[i].x<center.x && other[i].y>center.y)
		{
			left_down.push_back(other[i]);
		}
		//right_down
		if (other[i].x>center.x && other[i].y>center.y)
		{
			right_down.push_back(other[i]);
		}
	}
	refourvector.push_back(left_up);
	refourvector.push_back(left_down);
	refourvector.push_back(right_up);
	refourvector.push_back(right_down);
	return refourvector;
}

//find mini Point
int find_mini_Point(vector<cv::Point> Pointes, cv::Point center)
{
	if (Pointes.size() == 0)
	{
		cout << "没有点可以找了" << endl;
		return -1;
	}
	int reindex = 0;
	double length_mini = sqrt(pow(Pointes[0].x - center.x, 2) + pow(Pointes[0].y - center.y, 2));
	for (int i = 0; i<Pointes.size(); i++)
	{
		double length = sqrt(pow(Pointes[i].x - center.x, 2) + pow(Pointes[i].y - center.y, 2));
		if (length<length_mini)
		{
			length_mini = length;
			reindex = i;
		}
	}
	return reindex;
}

// 记录点的位置信息 删除已找到的点
cv::Point locationAnddelete(vector<cv::Point> & right_down, int index, vector<vector<double>> &right_down_location, int row, int col, cv::Point center,
	vector<double> row_margin, vector<double> col_margin, int mode /*mode=1,水平边查找；mode=2，竖直边查找，mode=3,斜边查找*/, int Quadrant/*象限*/, int threshold)
{
	vector<double> location_temp;
	cv::Point rePoint;
	cv::Point miss_point;
	//判断是否有缺点
	double length = sqrt(pow(right_down[index].x - center.x, 2) + pow(right_down[index].y - center.y, 2));
	double length_modle;
	int row_temp = row;
	int col_temp = col;
	if (mode == 1)
	{
		if (col > row_margin.size())
		{
			col_temp = row_margin.size();
		}
		length_modle = row_margin[row_margin.size() - col_temp];//水平边
	}
	if (mode == 2)
	{
		if (row > col_margin.size())
		{
			row_temp = col_margin.size();
		}
		length_modle = col_margin[col_margin.size() - row_temp];//竖直边
	}
	if (mode == 3)
	{
		if (row > col_margin.size())
		{
			row_temp = col_margin.size();
		}
		if (col > row_margin.size())
		{
			col_temp = row_margin.size();
		}
		length_modle = sqrt(pow(row_margin[row_margin.size() - col_temp], 2) + pow(col_margin[col_margin.size() - row_temp], 2));//斜边
	}

	if (abs(length - length_modle) > threshold)//说明点没找对，right_domn有缺少的点
	{
		cout << "缺少点了！！！" << endl;
		if (mode == 1)//水平边
		{
			//判断是哪一个象限
			if (Quadrant == 1)
			{
				miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y);
			}
			if (Quadrant == 2)
			{
				miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y);
			}
			if (Quadrant == 3)
			{
				miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y);
			}
			if (Quadrant == 4)
			{
				miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y);
			}

		}
		if (mode == 2)//竖直边
		{
			//判断是哪一个象限
			if (Quadrant == 1)
			{
				miss_point = cv::Point(center.x, center.y - col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 2)
			{
				miss_point = cv::Point(center.x, center.y - col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 3)
			{
				miss_point = cv::Point(center.x, center.y + col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 4)
			{
				miss_point = cv::Point(center.x, center.y + col_margin[col_margin.size() - row_temp]);
			}

		}
		if (mode == 3)//斜边
		{
			if (Quadrant == 1)
			{
				miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y - col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 2)
			{
				miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y - col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 3)
			{
				miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y + col_margin[col_margin.size() - row_temp]);
			}
			if (Quadrant == 4)
			{
				miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y + col_margin[col_margin.size() - row_temp]);
			}

		}
		//right_down_location 接收
		location_temp.push_back(miss_point.x);
		location_temp.push_back(miss_point.y);
		location_temp.push_back(row);
		location_temp.push_back(col);
		rePoint = miss_point;
	}
	else // right_down没有缺少的点
	{
		if (mode == 1)
		{
			//没有缺少点，但是没有找对（比如水平方向上缺少点的时候）
			if (abs(right_down[index].y - center.y) > 20)//水平方向上找点出现错误
			{
				cout << "水平方向找的时候，y值相差太大，没有找对！！！" << endl;
				//补充一个水平
				//判断是哪一个象限
				if (Quadrant == 1)
				{
					miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y);
				}
				if (Quadrant == 2)
				{
					miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y);
				}
				if (Quadrant == 3)
				{
					miss_point = cv::Point(center.x - row_margin[row_margin.size() - col_temp], center.y);
				}
				if (Quadrant == 4)
				{
					miss_point = cv::Point(center.x + row_margin[row_margin.size() - col_temp], center.y);
				}

				//right_down_location 接收
				location_temp.push_back(miss_point.x);
				location_temp.push_back(miss_point.y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = miss_point;
			}
			else
			{
				//right_down_location 接收
				location_temp.push_back(right_down[index].x);
				location_temp.push_back(right_down[index].y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = right_down[index];
				//删除找到的点
				right_down.erase(right_down.begin() + index);
			}

		}

		if (mode == 2)
		{
			if (abs(right_down[index].x - center.x) > 20)//水平方向上找点出现错误
			{
				cout << "竖直方向找的时候，x值相差太大，没有找对！！！" << endl;
				//判断是哪一个象限
				if (Quadrant == 1)
				{
					miss_point = cv::Point(center.x, center.y - col_margin[col_margin.size() - row_temp]);
				}
				if (Quadrant == 2)
				{
					miss_point = cv::Point(center.x, center.y - col_margin[col_margin.size() - row_temp]);
				}
				if (Quadrant == 3)
				{
					miss_point = cv::Point(center.x, center.y + col_margin[col_margin.size() - row_temp]);
				}
				if (Quadrant == 4)
				{
					miss_point = cv::Point(center.x, center.y + col_margin[col_margin.size() - row_temp]);
				}
				//right_down_location 接收
				location_temp.push_back(miss_point.x);
				location_temp.push_back(miss_point.y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = miss_point;
			}
			else
			{
				//right_down_location 接收
				location_temp.push_back(right_down[index].x);
				location_temp.push_back(right_down[index].y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = right_down[index];
				//删除找到的点
				right_down.erase(right_down.begin() + index);
			}

		}

		if (mode == 3)
		{
			//right_down_location 接收
			location_temp.push_back(right_down[index].x);
			location_temp.push_back(right_down[index].y);
			location_temp.push_back(row);
			location_temp.push_back(col);
			rePoint = right_down[index];
			//删除找到的点
			right_down.erase(right_down.begin() + index);
		}
	}
	right_down_location.push_back(location_temp);
	return rePoint;
}

int location(vector<vector<double>> &right_down_location, vector<cv::Point>right_down, cv::Point cameracenter_JZ, vector<double>RotaterowMargin_JZ, vector<double>RotatecolMargin_JZ,
	int Quadrant, int threshold)
{
	vector<cv::Point> row_Point;
	vector<cv::Point> col_Point;
	cv::Point home;
	cv::Point temp;


	//找到right_down_home点
	//vector<vector<double>> right_down_location;
	int row = 1;
	int col = 1;
	//找到 home点
	cv::Point right_down_home;
	int index = find_mini_Point(right_down, cameracenter_JZ);

	if (index == -1)
	{
		cout << "起始点home没有找到！！！" << endl;
		return -2;
	}
	right_down_home = right_down[index];
	locationAnddelete(right_down, index, right_down_location, row, col, cameracenter_JZ, RotaterowMargin_JZ, RotatecolMargin_JZ, 3, Quadrant, threshold);

	int i = 0;
	while (right_down.size() != 0)
	{
		i += 1;
		//第一次需要从right_down_home点开始找三个点（row，col，home）；
		if (i == 1)
		{
			index = find_mini_Point(right_down, right_down_home);
			if (index == -1)
			{
				cout << "第一次找最小点，没有找到" << endl;
				break;
			}
			cv::Point temp_one = right_down[index];
			if (abs(temp_one.x - right_down_home.x)<20) //说明点在上边或者下边
			{
				temp = locationAnddelete(right_down, index, right_down_location, row + 1, col, right_down_home, RotaterowMargin_JZ, RotatecolMargin_JZ, 2, Quadrant, threshold);
				row_Point.push_back(temp);

				index = find_mini_Point(right_down, right_down_home);
				if (index == -1)
				{
					cout << "第一次找最小点，没有找到" << endl;
					break;
				}
				temp = locationAnddelete(right_down, index, right_down_location, row, col + 1, right_down_home, RotaterowMargin_JZ, RotatecolMargin_JZ, 1, Quadrant, threshold);
				col_Point.push_back(temp);
			}
			else // 在右边
			{
				temp = locationAnddelete(right_down, index, right_down_location, row, col + 1, right_down_home, RotaterowMargin_JZ, RotatecolMargin_JZ, 1, Quadrant, threshold);
				col_Point.push_back(temp);

				index = find_mini_Point(right_down, right_down_home);
				if (index == -1)
				{
					cout << "第一次找最小点，没有找到" << endl;
					break;
				}
				temp = locationAnddelete(right_down, index, right_down_location, row + 1, col, right_down_home, RotaterowMargin_JZ, RotatecolMargin_JZ, 2, Quadrant, threshold);
				row_Point.push_back(temp);
			}
			//查找home
			index = find_mini_Point(right_down, right_down_home);
			if (index == -1)
			{
				cout << "第一次找最小点，没有找到" << endl;
				break;
			}
			temp = locationAnddelete(right_down, index, right_down_location, row + 1, col + 1, right_down_home, RotaterowMargin_JZ, RotatecolMargin_JZ, 3, Quadrant, threshold);
			home = temp;
			//每一次找完之后，更新(x,y),便于下一次使用
			row += 1;
			col += 1;
		}


		vector<cv::Point> temp_row_Point;
		vector<cv::Point> temp_col_Point;
		//其他
		//1. 通过row_Point 找下一个水平上的点
		if (row_Point.size() == 0)
		{
			cout << "row_Point is null" << endl;
			break;
		}
		for (int k = 0; k<row_Point.size(); k++)
		{
			if (k == 0)//第一个点
			{
				//第一次
				index = find_mini_Point(right_down, row_Point[0]);

				if (index == -1)
				{
					cout << "row的第一个点，没有找到最小点" << endl;
					break;
				}
				temp = locationAnddelete(right_down, index, right_down_location, row + 1, k + 1, row_Point[0], RotaterowMargin_JZ, RotatecolMargin_JZ, 2, Quadrant, threshold);
				temp_row_Point.push_back(temp);
			}
			//其他水平方向上的点
			index = find_mini_Point(right_down, row_Point[k]);
			if (index == -1)
			{
				cout << "row---没有找到最小点" << endl;
				break;
			}
			temp = locationAnddelete(right_down, index, right_down_location, row + 1, k + 2, row_Point[k], RotaterowMargin_JZ, RotatecolMargin_JZ, 3, Quadrant, threshold);
			temp_row_Point.push_back(temp);
		}

		//2. 通过col_Point 找下一个竖直方向上的点
		if (col_Point.size() == 0)
		{
			cout << "col_Point is null" << endl;
			break;
		}
		for (int j = 0; j<col_Point.size(); j++)
		{
			if (j == 0)//第一个点
			{
				index = find_mini_Point(right_down, col_Point[0]);
				if (index == -1)
				{
					cout << "col第一个点，没有找到最小点" << endl;
					break;
				}
				temp = locationAnddelete(right_down, index, right_down_location, j + 1, col + 1, col_Point[0], RotaterowMargin_JZ, RotatecolMargin_JZ, 1, Quadrant, threshold);
				temp_col_Point.push_back(temp);
			}
			index = find_mini_Point(right_down, col_Point[j]);
			if (index == -1)
			{
				cout << "col---没有找到最小点" << endl;
				break;
			}
			temp = locationAnddelete(right_down, index, right_down_location, j + 2, col + 1, col_Point[j], RotaterowMargin_JZ, RotatecolMargin_JZ, 3, Quadrant, threshold);
			temp_col_Point.push_back(temp);
		}
		//3. 通过home 斜向推进

		index = find_mini_Point(right_down, home);
		if (index == -1)
		{
			cout << "home---没有找到最小点" << endl;
			break;
		}
		temp = locationAnddelete(right_down, index, right_down_location, row + 1, col + 1, home, RotaterowMargin_JZ, RotatecolMargin_JZ, 3, Quadrant, threshold);
		home = temp;
		row += 1;
		col += 1;
		row_Point.clear();
		row_Point = temp_row_Point;
		col_Point.clear();
		col_Point = temp_col_Point;

	}
	return 100;
}

//配置ori位置信息
vector<vector<double>> location_ori(vector<cv::Point> Fit_othercooder, vector<vector<double>> fit_location, vector<cv::Point> othercooder, int&rownum, int &colnum)
{
	vector<vector<double>> re_location;
	for (int i = 0; i <Fit_othercooder.size(); i++)
	{
		for (int j = 0; j < fit_location.size(); j++)
		{
			if (Fit_othercooder[i].x == fit_location[j][0] &&
				Fit_othercooder[i].y == fit_location[j][1])
			{
				vector<double> temp;
				temp.push_back(othercooder[i].x);
				temp.push_back(othercooder[i].y);
				temp.push_back(fit_location[j][2]);
				temp.push_back(fit_location[j][3]);
				// 为了得到最初拟合图像中又多少列多上行，也就是原始图像中的row，col
				if (fit_location[j][3] > rownum)
				{
					rownum = fit_location[j][3];
				}
				if (fit_location[j][2] > colnum)
				{
					colnum = fit_location[j][2];
				}
				re_location.push_back(temp);
				break;
			}
		}
	}
	return re_location;
}

// 生成完美的眼睛图像
vector<vector<double>> generate_Eye_other(vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, int Quadrant, int rownum, int colnum)
{
	vector<vector<double>> location_Eye;
	cv::Point col_home;
	cv::Point row_home;
	cv::Point home = center;


	for (int i = 0; i < colnum; i++)
	{
		int y_index = R_cameracol_JZ.size() - i - 1;
		if (y_index < 0)
		{
			y_index = 0;
		}
		if (Quadrant == 1)
		{
			col_home.y = home.y - R_cameracol_JZ[y_index];
			col_home.x = home.x;
			home.y = col_home.y;
		}
		if (Quadrant == 2)
		{
			col_home.y = home.y - R_cameracol_JZ[y_index];
			col_home.x = home.x;
			home.y = col_home.y;
		}
		if (Quadrant == 3)
		{
			col_home.y = home.y + R_cameracol_JZ[y_index];
			col_home.x = home.x;
			home.y = col_home.y;
		}
		if (Quadrant == 4)
		{
			col_home.y = home.y + R_cameracol_JZ[y_index];
			col_home.x = home.x;
			home.y = col_home.y;
		}
		for (int j = 0; j < rownum; j++)
		{
			vector<double> temp_vector;
			int x_index = R_camerarow_JZ.size() - j - 1;
			if (x_index < 0)
			{
				x_index = 0;
			}
			if (Quadrant == 1)
			{
				row_home.x = col_home.x + R_camerarow_JZ[x_index];
				col_home.x = row_home.x;
				row_home.y = col_home.y;

			}
			if (Quadrant == 2)
			{
				row_home.x = col_home.x - R_camerarow_JZ[x_index];
				col_home.x = row_home.x;
				row_home.y = col_home.y;
			}
			if (Quadrant == 3)
			{
				row_home.x = col_home.x - R_camerarow_JZ[x_index];
				col_home.x = row_home.x;
				row_home.y = col_home.y;
			}
			if (Quadrant == 4)
			{
				row_home.x = col_home.x + R_camerarow_JZ[x_index];
				col_home.x = row_home.x;
				row_home.y = col_home.y;
			}
			temp_vector.push_back(row_home.x);
			temp_vector.push_back(row_home.y);
			temp_vector.push_back(i + 1);
			temp_vector.push_back(j + 1);
			location_Eye.push_back(temp_vector);
		}
	}
	return location_Eye;
}
//生成eye图像的row col
void generate_row_col(int colnum, int rownum, vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, vector<cv::Point> &rowcooder_eye, vector<cv::Point>& colcooder_eye)
{
	cv::Point home1 = center, home2 = center;
	for (int i = 0; i < colnum; i++)
	{
		int y_index = R_cameracol_JZ.size() - i - 1;
		if (y_index < 0)
		{
			y_index = 0;
		}
		cv::Point temp1, temp2;
		temp1.x = center.x;
		temp1.y = home1.y - R_cameracol_JZ[y_index];
		home1.y = temp1.y;

		temp2.x = center.x;
		temp2.y = home2.y + R_cameracol_JZ[y_index];
		home2.y = temp2.y;
		colcooder_eye.push_back(temp1);
		colcooder_eye.push_back(temp2);
	}

	home1 = center, home2 = center;
	for (int i = 0; i < rownum; i++)
	{
		int x_index = R_camerarow_JZ.size() - i - 1;
		if (x_index < 0)
		{
			x_index = 0;
		}
		cv::Point temp1;
		temp1.y = center.y;
		temp1.x = home1.x - R_camerarow_JZ[x_index];
		home1.x = temp1.x;

		cv::Point temp2;
		temp2.y = center.y;
		temp2.x = home2.x + R_camerarow_JZ[x_index];
		home2.x = temp2.x;

		rowcooder_eye.push_back(temp1);
		rowcooder_eye.push_back(temp2);
	}
	colcooder_eye.push_back(center);
	rowcooder_eye.push_back(center);

}

//根据camimage图像的坐标调整eyeimage的坐标
void moveOri(vector<vector<double>> &right_down_location_ori, vector<vector<double>> right_down_location, vector<vector<double>> right_down_location_Eye)
{
	int cam_index = -1;
	int eye_index = -1;
	for (int i = 0; i < right_down_location_ori.size(); i++)
	{
		for (int j = 0; j < right_down_location.size(); j++)
		{
			if (right_down_location_ori[i][2] == right_down_location[j][2] &&
				right_down_location_ori[i][3] == right_down_location[j][3])
			{
				cam_index = j;
				break;
			}
		}
		for (int k = 0; k < right_down_location_Eye.size(); k++)
		{
			if (right_down_location_ori[i][2] == right_down_location_Eye[k][2] &&
				right_down_location_ori[i][3] == right_down_location_Eye[k][3])
			{
				eye_index = k;
				break;
			}
		}
		if (eye_index == -1 || cam_index == -1)
		{
			cout << "没找到对应的点" << endl;
			break;
		}
		//判断camImage上的点与eyeimage上的点到底相差多少
		double x_length = right_down_location[cam_index][0] - right_down_location_Eye[eye_index][0];
		double y_length = right_down_location[cam_index][1] - right_down_location_Eye[eye_index][1];
		if (sqrt(pow(x_length, 2) + pow(y_length, 2))<5)
		{
			cout << "这个点不用调整！！！" << endl;
		}
		else
		{
			cout << "参与调整了" << endl;
			
			if (x_length<0)
			{
				right_down_location_ori[i][0] += 1;
			}
			else
			{
				right_down_location_ori[i][0] -= 1;
			}

			if (y_length<0)
			{
				right_down_location_ori[i][1] += 1;
			}
			else
			{
				right_down_location_ori[i][1] -= 1;
			}
			
		}
	}
}



void ori_fit(cv::Mat oriImage, vector<cv::Point> &othercooder, vector<cv::Point> &rowcooder, vector<cv::Point> &colcooder, cv::Point &center, vector<cv::Point>& Fit_othercooder,
	vector<vector<double>> &right_down_location_fit, vector<vector<double>> &left_up_location_fit, vector<vector<double>> &left_down_location_fit,vector<vector<double>> &right_up_location_fit)
{
	int font_face = cv::FONT_HERSHEY_COMPLEX;

	cv::Mat oriChannel[3];
	cv::split(oriImage, oriChannel);


	/************************************ Ori_Image ***********************************************/
	// 1. red channel
	cv::Mat oriRed;
	oriChannel[2].copyTo(oriRed);
	// get colcooder
	
	getcoorder(oriRed, colcooder);
	colsort(colcooder, 1);

	// 2. bule channel
	cv::Mat oriBlue;
	oriChannel[0].copyTo(oriBlue);
	// get rowcooder
	
	getcoorder(oriBlue, rowcooder);
	colsort(rowcooder, 0);


	// get center
	center  = getcenter(oriRed);

	//4. green channel
	cv::Mat oriGreen;
	oriChannel[1].copyTo(oriGreen);
	// get othercooder
	
	getcoorder(oriGreen, othercooder);


	//5. get margin of row and col
	vector<double> rowMargin;
	rowMargin = getMargin(rowcooder, center, 0);
	vector<double> colMargin;
	colMargin = getMargin(colcooder, center, 1);



	/************************************ Fit_Image ***********************************************/


	cv::Point fitCenter = cv::Point(900, 1000);
	//1. get fit_row cooder
	vector<cv::Point> Fit_rowcooder;
	Fit_rowcooder = getFitData(rowcooder, center, fitCenter);//

	//2. get fit_col cooder
	vector<cv::Point> Fit_colcooder;
	Fit_colcooder = getFitData(colcooder, center, fitCenter);//

	//3. get fit_other cooder
	
	Fit_othercooder = getFitData(othercooder, center, fitCenter);//

	//4. get fit_row_margin
	vector<double> Fit_rowMargin;
	Fit_rowMargin = getMargin(Fit_rowcooder, fitCenter, 0);

	//5. get col_row_margin
	vector<double> Fit_colMargin;
	Fit_colMargin = getMargin(Fit_colcooder, fitCenter, 1);

	// write fit_image 
	cv::Mat fitImage(1800, 2000, CV_8UC3, cv::Scalar::all(0));

	writefitImage(fitImage, Fit_rowcooder, 1);
	writefitImage(fitImage, Fit_colcooder, 0);
	writefitImage(fitImage, Fit_othercooder, 2);
	cv::circle(fitImage, fitCenter, 15, cv::Scalar(255, 0, 255), -1, 8);
	cv::imwrite("fitimage.jpg", fitImage);


	////////////定位拟合图像/////////////
	//Fit_othercooder 分成四个部分进行定位
	vector<vector<cv::Point>> refourvector_fit;
	refourvector_fit = separatefour(Fit_othercooder, fitCenter);

	vector<cv::Point> left_up_fit, left_domn_fit, right_up_fit, right_down_fit;
	left_up_fit = refourvector_fit[0];
	left_domn_fit = refourvector_fit[1];
	right_up_fit = refourvector_fit[2];
	right_down_fit = refourvector_fit[3];


	
	location(right_up_location_fit, right_up_fit, fitCenter, Fit_rowMargin, Fit_colMargin, 1, 20);//第一象限
	location(left_up_location_fit, left_up_fit, fitCenter, Fit_rowMargin, Fit_colMargin, 2, 20);//第二象限
	location(left_down_location_fit, left_domn_fit, fitCenter, Fit_rowMargin, Fit_colMargin, 3, 20);//第三象限
	location(right_down_location_fit, right_down_fit, fitCenter, Fit_rowMargin, Fit_colMargin, 4, 20);//第四象限


	//写位置信息（坐标编码）
	for (int i = 0; i < right_down_location_fit.size(); i++)
	{
		string text = doubleToString(right_down_location_fit[i][2]).append(",").append(doubleToString(right_down_location_fit[i][3]));
		cv::putText(fitImage, text, cv::Point(right_down_location_fit[i][0] - 5, right_down_location_fit[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_up_location_fit.size(); i++)
	{
		string text = doubleToString(left_up_location_fit[i][2]).append(",").append(doubleToString(left_up_location_fit[i][3]));
		cv::putText(fitImage, text, cv::Point(left_up_location_fit[i][0] - 5, left_up_location_fit[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_down_location_fit.size(); i++)
	{
		string text = doubleToString(left_down_location_fit[i][2]).append(",").append(doubleToString(left_down_location_fit[i][3]));
		cv::putText(fitImage, text, cv::Point(left_down_location_fit[i][0] - 5, left_down_location_fit[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < right_up_location_fit.size(); i++)
	{
		string text = doubleToString(right_up_location_fit[i][2]).append(",").append(doubleToString(right_up_location_fit[i][3]));
		cv::putText(fitImage, text, cv::Point(right_up_location_fit[i][0] - 5, right_up_location_fit[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	cv::imwrite("fitimage.jpg", fitImage);
}


void camera_location(cv::Mat cameraImage, vector<double> &RotatecolMargin_JZ, vector<double> &RotaterowMargin_JZ, cv::Point &cameracenter_JZ,
	vector<vector<double>>& right_down_location, vector<vector<double>>& left_up_location, vector<vector<double>>& left_down_location, vector<vector<double>>& right_up_location)
{
	int font_face = cv::FONT_HERSHEY_COMPLEX;

	/************************************ Camera_Image ***********************************************/
	//1.read image
	//cv::Mat cameraImage;
	//cameraImage = cv::imread("5.bmp");//
	//image_Camera.copyTo(cameraImage); //获取相机图像


	if (cameraImage.data == NULL)
	{
		cout << "cameraimage is null" << endl;
		
	}
	cv::imwrite("cameraImage.bmp", cameraImage);

	//2. image process
	ContrastImage(cameraImage, 2, -100);
	cv::Mat cameraChannel[3];
	cv::split(cameraImage, cameraChannel);
	static cv::Mat element1 = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(15, 15));
	static cv::Mat element2 = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));

	//3. red channel
	cv::Mat cameraRed;
	cameraChannel[2].copyTo(cameraRed);

	cv::threshold(cameraRed, cameraRed, 100, 255, 0);
	//cv::imwrite("test.jpg", cameraRed);
	cv::erode(cameraRed, cameraRed, element1);
	cv::dilate(cameraRed, cameraRed, element2);
	//get camera red channel
	vector<cv::Point> camera_Redcooder;
	getcoorder(cameraRed, camera_Redcooder);
	cv::Point cameraCenter = getcenter(cameraRed);//
	//cv::imwrite("test.jpg", cameraRed);

	//4. bule channel
	cv::Mat cameraBlue;
	cameraChannel[0].copyTo(cameraBlue);
	cv::threshold(cameraBlue, cameraBlue, 100, 255, 0);
	cv::erode(cameraBlue, cameraBlue, element1);
	cv::dilate(cameraBlue, cameraBlue, element2);
	//cv::imwrite("test.jpg", cameraBlue);
	// get bule channel
	vector<cv::Point> camera_Bluecooder;
	getcoorder(cameraBlue, camera_Bluecooder); //

	//5. green channel 
	cv::Mat cameraGreen;
	cameraChannel[1].copyTo(cameraGreen);
	cv::threshold(cameraGreen, cameraGreen, 10, 255, 0);
	//cv::imwrite("test.jpg", cameraGreen);

	//cv::erode(cameraGreen, cameraGreen, element1);
	//cv::dilate(cameraGreen, cameraGreen, element2);

	//cv::imwrite("test.jpg", cameraGreen);
	// get green channel
	vector<cv::Point> camera_Greencooder;
	getcoorder(cameraGreen, camera_Greencooder); //

	//6. all channel pointes need sure again!!!
	vector<cv::Point>camera_colcooder;
	camera_colcooder = cameraCooderRight(camera_Redcooder, cameraCenter, 1); //camera_col_pointes

	vector<cv::Point>camera_rowcooder;
	camera_rowcooder = cameraCooderRight(camera_Bluecooder, cameraCenter, 0);  //camera_row_pointes

	vector<cv::Point>camera_othercooder;
	camera_othercooder = otherCooderRight(camera_rowcooder, camera_colcooder, camera_Greencooder);//camer_other_pointes

	// write camra image
	cv::Mat prefectcameraImage(cameraImage.size(), CV_8UC3, cv::Scalar(0, 0, 0));
	writefitImage(prefectcameraImage, camera_rowcooder, 1);
	writefitImage(prefectcameraImage, camera_colcooder, 0);
	writefitImage(prefectcameraImage, camera_othercooder, 2);

	cv::circle(prefectcameraImage, cameraCenter, 15, cv::Scalar(255, 0, 255), -1, 8);
	cv::imwrite("prefectcameraImage.jpg", prefectcameraImage);



	/******************* write csv file ***************/




	/************************************ JZ Camera_Image ***********************************************/

	//distortion camera image 
	cv::Mat intrinsicMatrix = cv::Mat::zeros(3, 3, CV_64F);
	cv::Mat distortionCoeff = cv::Mat::ones(1, 4, CV_64F);

	intrinsicMatrix.at<double>(0, 0) = 1601.575552125639;
	intrinsicMatrix.at<double>(0, 1) = 0;
	intrinsicMatrix.at<double>(0, 2) = 2047.217820040776;
	intrinsicMatrix.at<double>(1, 0) = 0;
	intrinsicMatrix.at<double>(1, 1) = 1602.774750782801;
	intrinsicMatrix.at<double>(1, 2) = 1481.768466397291;
	intrinsicMatrix.at<double>(2, 0) = 0;
	intrinsicMatrix.at<double>(2, 1) = 0;
	intrinsicMatrix.at<double>(2, 2) = 1;

	distortionCoeff.at<double>(0, 0) = -0.00730083;
	distortionCoeff.at<double>(0, 1) = -0.000315583;
	distortionCoeff.at<double>(0, 2) = -0.000565603;
	distortionCoeff.at<double>(0, 3) = -0.000978041;


	vector<cv::Point> cameracenter_right;
	cameracenter_right.push_back(cameraCenter);

	// get JZ image center pointe
	vector<cv::Point> cameracenter_JZ_temp;
	myUndistortPoints(cameracenter_right, cameracenter_JZ_temp, intrinsicMatrix, distortionCoeff);
	cameracenter_JZ = cameracenter_JZ_temp[0];

	//get JZ image row cooder
	vector<cv::Point> camerarow_JZ;
	myUndistortPoints(camera_rowcooder, camerarow_JZ, intrinsicMatrix, distortionCoeff);

	//get JZ image col cooder
	vector<cv::Point> cameracol_JZ;
	myUndistortPoints(camera_colcooder, cameracol_JZ, intrinsicMatrix, distortionCoeff);

	// get JZ image other cooder
	vector<cv::Point> cameraother_JZ;
	myUndistortPoints(camera_othercooder, cameraother_JZ, intrinsicMatrix, distortionCoeff);

	// sort JZ row cooder and col cooder
	colsort(camerarow_JZ, 0);
	colsort(cameracol_JZ, 1);

	//讲cameracenter_JZ camerarow_JZ,cameracol_JZ,cameraother_JZ 全部平移
	for (int i = 0; i < camerarow_JZ.size(); i++)
	{
		camerarow_JZ[i].x += 2500 - cameracenter_JZ.x;
		camerarow_JZ[i].y += 2000 - cameracenter_JZ.y;
	}
	for (int i = 0; i < cameracol_JZ.size(); i++)
	{
		cameracol_JZ[i].x += 2500 - cameracenter_JZ.x;
		cameracol_JZ[i].y += 2000 - cameracenter_JZ.y;
	}
	for (int i = 0; i < cameraother_JZ.size(); i++)
	{
		cameraother_JZ[i].x += 2500 - cameracenter_JZ.x;
		cameraother_JZ[i].y += 2000 - cameracenter_JZ.y;
	}
	cameracenter_JZ.x += 2500 - cameracenter_JZ.x;
	cameracenter_JZ.y += 2000 - cameracenter_JZ.y;

	// write JZ image
	cv::Mat cameraImage_JZ(4000, 5000, CV_8UC3, cv::Scalar(0, 0, 0));
	writefitImage(cameraImage_JZ, camerarow_JZ, 1);
	writefitImage(cameraImage_JZ, cameracol_JZ, 0);
	writefitImage(cameraImage_JZ, cameraother_JZ, 2);
	cv::circle(cameraImage_JZ, cameracenter_JZ, 15, cv::Scalar(255, 0, 255), -1, 8);
	cv::imwrite("DESprefectcameraImage.jpg", cameraImage_JZ);



	/************************************ Rotate JZ_Camera_Image ***********************************************/

	// rotate JZ image
	double angle1, angle2, angle3, angle4, angle;

	double temp_left = atan2((camerarow_JZ[0].y - cameracenter_JZ.y), (cameracenter_JZ.x - camerarow_JZ[0].x));
	angle1 = temp_left;

	double temp_up = atan2((cameracol_JZ[0].x - cameracenter_JZ.x), (cameracenter_JZ.y - cameracol_JZ[0].y));
	angle2 = -temp_up;

	double temp_right = atan2((camerarow_JZ[camerarow_JZ.size() - 1].y - cameracenter_JZ.y), (camerarow_JZ[camerarow_JZ.size() - 1].x - cameracenter_JZ.x));
	angle3 = -temp_right;

	double temp_down = atan2((cameracol_JZ[cameracol_JZ.size() - 1].x - cameracenter_JZ.x), (cameracol_JZ[cameracol_JZ.size() - 1].y - cameracenter_JZ.y));
	angle4 = temp_down;

	angle = -(angle1 + angle2 + angle3 + angle4) / 4;

	vector<cv::Point> R_camerarow_JZ;
	R_camerarow_JZ = CalrotatePoint(camerarow_JZ, cameracenter_JZ, angle);

	vector<cv::Point> R_cameracol_JZ;
	R_cameracol_JZ = CalrotatePoint(cameracol_JZ, cameracenter_JZ, angle);

	vector<cv::Point> R_cameraother_JZ;
	R_cameraother_JZ = CalrotatePoint(cameraother_JZ, cameracenter_JZ, angle);


	
	RotaterowMargin_JZ = getMargin(R_camerarow_JZ, cameracenter_JZ, 0);
	
	RotatecolMargin_JZ = getMargin(R_cameracol_JZ, cameracenter_JZ, 1);


	cv::Mat R_cameraImage_JZ(4000, 5000, CV_8UC3, cv::Scalar(0, 0, 0));
	writefitImage(R_cameraImage_JZ, R_camerarow_JZ, 1);
	writefitImage(R_cameraImage_JZ, R_cameracol_JZ, 0);
	writefitImage(R_cameraImage_JZ, R_cameraother_JZ, 2);
	cv::circle(R_cameraImage_JZ, cameracenter_JZ, 15, cv::Scalar(255, 0, 255), -1, 8);
	cv::imwrite("R_DESprefectcameraImage.jpg", R_cameraImage_JZ);


	/////////////////定位旋转图像//////////////////
	//Rotatecameraother_JZ 旋转之后的other点，中心点：cameracenter_JZ
	vector<vector<cv::Point>> refourvector;
	refourvector = separatefour(R_cameraother_JZ, cameracenter_JZ);

	vector<cv::Point> left_up, left_domn, right_up, right_down;
	left_up = refourvector[0];
	left_domn = refourvector[1];
	right_up = refourvector[2];
	right_down = refourvector[3];



	location(right_up_location, right_up, cameracenter_JZ, RotaterowMargin_JZ, RotatecolMargin_JZ, 1, 61);//第一象限
	location(left_up_location, left_up, cameracenter_JZ, RotaterowMargin_JZ, RotatecolMargin_JZ, 2, 61);//第二象限
	location(left_down_location, left_domn, cameracenter_JZ, RotaterowMargin_JZ, RotatecolMargin_JZ, 3, 61);//第三象限
	location(right_down_location, right_down, cameracenter_JZ, RotaterowMargin_JZ, RotatecolMargin_JZ, 4, 61);//第四象限

	for (int i = 0; i < right_down_location.size(); i++)
	{
		string text = doubleToString(right_down_location[i][2]).append(",").append(doubleToString(right_down_location[i][3]));
		cv::putText(R_cameraImage_JZ, text, cv::Point(right_down_location[i][0] - 5, right_down_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_up_location.size(); i++)
	{
		string text = doubleToString(left_up_location[i][2]).append(",").append(doubleToString(left_up_location[i][3]));
		cv::putText(R_cameraImage_JZ, text, cv::Point(left_up_location[i][0] - 5, left_up_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_down_location.size(); i++)
	{
		string text = doubleToString(left_down_location[i][2]).append(",").append(doubleToString(left_down_location[i][3]));
		cv::putText(R_cameraImage_JZ, text, cv::Point(left_down_location[i][0] - 5, left_down_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < right_up_location.size(); i++)
	{
		string text = doubleToString(right_up_location[i][2]).append(",").append(doubleToString(right_up_location[i][3]));
		cv::putText(R_cameraImage_JZ, text, cv::Point(right_up_location[i][0] - 5, right_up_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	cv::imwrite("R_DESprefectcameraImage.jpg", R_cameraImage_JZ);

}



void corner_cor(vector<cv::Point> &Fit_othercooder, vector<cv::Point> &othercooder, vector<double> &RotatecolMargin_JZ, vector<double> &RotaterowMargin_JZ, cv::Point &cameracenter_JZ,
	vector<vector<double>> &right_down_location_fit, vector<vector<double>> &left_up_location_fit, vector<vector<double>> &left_down_location_fit, vector<vector<double>> &right_up_location_fit,
	vector<vector<double>> &right_down_location_ori, vector<vector<double>>& left_up_location_ori, vector<vector<double>>& left_down_location_ori, vector<vector<double>>& right_up_location_ori,
	vector<vector<double>> & right_down_location_Eye, vector<vector<double>>& left_up_location_Eye, vector<vector<double>>& left_down_location_Eye, vector<vector<double>>& right_up_location_Eye)
{
	/**************************************************** 角点对应 *****************************************************/
	//////////定位ori图像//////////////
	int font_face = cv::FONT_HERSHEY_COMPLEX;

	int ori_rownum = 0;
	int ori_colnum = 0;
	right_down_location_ori = location_ori(Fit_othercooder, right_down_location_fit, othercooder, ori_rownum, ori_colnum);
	left_up_location_ori = location_ori(Fit_othercooder, left_up_location_fit, othercooder, ori_rownum, ori_colnum);
	left_down_location_ori = location_ori(Fit_othercooder, left_down_location_fit, othercooder, ori_rownum, ori_colnum);
	right_up_location_ori = location_ori(Fit_othercooder, right_up_location_fit, othercooder, ori_rownum, ori_colnum);


	/************************************************* Generate prefect Eye_Image ***********************************************/


	//生成需要“横平竖直“的图像坐标及位置信息（在R_cameraImage_JZ的基础上）	
	// R_cameraImage_JZ: right_down_location, left_up_location, left_down_location, right_up_location

	
	cv::Mat prefectEyeImage(4000, 5000, CV_8UC3, cv::Scalar::all(0));

	//生成prefecteye 图像的参数；
	vector<cv::Point> rowcooder_eye, colcooder_eye;
	generate_row_col(ori_colnum, ori_rownum, RotaterowMargin_JZ, RotatecolMargin_JZ, cameracenter_JZ, rowcooder_eye, colcooder_eye);
	writefitImage(prefectEyeImage, rowcooder_eye, 1);
	writefitImage(prefectEyeImage, colcooder_eye, 0);
	cv::Point center_eye = cameracenter_JZ;

	//四部分点的位置信息
	right_up_location_Eye = generate_Eye_other(RotaterowMargin_JZ, RotatecolMargin_JZ, cameracenter_JZ, 1, ori_rownum, ori_colnum);
	left_up_location_Eye = generate_Eye_other(RotaterowMargin_JZ, RotatecolMargin_JZ, cameracenter_JZ, 2, ori_rownum, ori_colnum);
	left_down_location_Eye = generate_Eye_other(RotaterowMargin_JZ, RotatecolMargin_JZ, cameracenter_JZ, 3, ori_rownum, ori_colnum);
	right_down_location_Eye = generate_Eye_other(RotaterowMargin_JZ, RotatecolMargin_JZ, cameracenter_JZ, 4, ori_rownum, ori_colnum);

	for (int i = 0; i < right_up_location_Eye.size(); i++)
	{
		cv::circle(prefectEyeImage, cv::Point(right_up_location_Eye[i][0], right_up_location_Eye[i][1]), 5, cv::Scalar(0, 255, 0), -1, 8);
		string text = doubleToString(right_up_location_Eye[i][2]).append(",").append(doubleToString(right_up_location_Eye[i][3]));
		cv::putText(prefectEyeImage, text, cv::Point(right_up_location_Eye[i][0] - 5, right_up_location_Eye[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_up_location_Eye.size(); i++)
	{
		cv::circle(prefectEyeImage, cv::Point(left_up_location_Eye[i][0], left_up_location_Eye[i][1]), 5, cv::Scalar(0, 255, 0), -1, 8);
		string text = doubleToString(left_up_location_Eye[i][2]).append(",").append(doubleToString(left_up_location_Eye[i][3]));
		cv::putText(prefectEyeImage, text, cv::Point(left_up_location_Eye[i][0] - 5, left_up_location_Eye[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_down_location_Eye.size(); i++)
	{
		cv::circle(prefectEyeImage, cv::Point(left_down_location_Eye[i][0], left_down_location_Eye[i][1]), 5, cv::Scalar(0, 255, 0), -1, 8);
		string text = doubleToString(left_down_location_Eye[i][2]).append(",").append(doubleToString(left_down_location_Eye[i][3]));
		cv::putText(prefectEyeImage, text, cv::Point(left_down_location_Eye[i][0] - 5, left_down_location_Eye[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < right_down_location_Eye.size(); i++)
	{
		cv::circle(prefectEyeImage, cv::Point(right_down_location_Eye[i][0], right_down_location_Eye[i][1]), 5, cv::Scalar(0, 255, 0), -1, 8);
		string text = doubleToString(right_down_location_Eye[i][2]).append(",").append(doubleToString(right_down_location_Eye[i][3]));
		cv::putText(prefectEyeImage, text, cv::Point(right_down_location_Eye[i][0] - 5, right_down_location_Eye[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	cv::circle(prefectEyeImage, center_eye, 15, cv::Scalar(255, 0, 255), -1, 8);
	cv::imwrite("prefectEyeImage.jpg", prefectEyeImage);

	
}

void write_ori_location(vector<vector<double>> &right_down_location_ori, vector<vector<double>>& left_up_location_ori, vector<vector<double>>& left_down_location_ori, vector<vector<double>>& right_up_location_ori)
{
	int font_face = cv::FONT_HERSHEY_COMPLEX;
	cv::Mat Ori_location = cv::imread("Ori_location.bmp");
	for (int i = 0; i < right_down_location_ori.size(); i++)
	{
		string text = doubleToString(right_down_location_ori[i][2]).append(",").append(doubleToString(right_down_location_ori[i][3]));
		cv::putText(Ori_location, text, cv::Point(right_down_location_ori[i][0] - 5, right_down_location_ori[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_up_location_ori.size(); i++)
	{
		string text = doubleToString(left_up_location_ori[i][2]).append(",").append(doubleToString(left_up_location_ori[i][3]));
		cv::putText(Ori_location, text, cv::Point(left_up_location_ori[i][0] - 5, left_up_location_ori[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < left_down_location_ori.size(); i++)
	{
		string text = doubleToString(left_down_location_ori[i][2]).append(",").append(doubleToString(left_down_location_ori[i][3]));
		cv::putText(Ori_location, text, cv::Point(left_down_location_ori[i][0] - 5, left_down_location_ori[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	for (int i = 0; i < right_up_location_ori.size(); i++)
	{
		string text = doubleToString(right_up_location_ori[i][2]).append(",").append(doubleToString(right_up_location_ori[i][3]));
		cv::putText(Ori_location, text, cv::Point(right_up_location_ori[i][0] - 5, right_up_location_ori[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
	}
	cv::imwrite("Ori_location.bmp", Ori_location);
}

void reori_show(cv::Mat oriImage, vector<cv::Point> &rowcooder, vector<cv::Point> &colcooder, cv::Point &center,
	vector<vector<double>> &right_down_location_ori, vector<vector<double>> &left_up_location_ori, vector<vector<double>> &left_down_location_ori, vector<vector<double>> &right_up_location_ori,
	vector<vector<double>>  &right_down_location_Eye, vector<vector<double>> &left_up_location_Eye, vector<vector<double>> &left_down_location_Eye, vector<vector<double>> &right_up_location_Eye,
	vector<vector<double>> &right_down_location, vector<vector<double>> &left_up_location, vector<vector<double>> &left_down_location, vector<vector<double>> &right_up_location)
{
	/*
	R_cameraImage_JZ:
	right_down_location, left_up_location, left_down_location, right_up_location
	eyeImage:
	right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye;
	oriImage:
	right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori;
	*/
	//调整原始图片
	moveOri(right_down_location_ori, right_down_location, right_down_location_Eye);
	moveOri(left_up_location_ori, left_up_location, left_up_location_Eye);
	moveOri(left_down_location_ori, left_down_location, left_down_location_Eye);
	moveOri(right_up_location_ori, right_up_location, right_up_location_Eye);
	//重新生成原始图片
	cv::Mat reOriImage(oriImage.size(), CV_8UC3, cv::Scalar::all(0));
	writefitImage(reOriImage, rowcooder, 1);
	writefitImage(reOriImage, colcooder, 0);
	cv::circle(reOriImage, center, 15, cv::Scalar(255, 0, 255), -1, 8);

	for (int i = 0; i < right_down_location_ori.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(right_down_location_ori[i][0], right_down_location_ori[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}

	for (int i = 0; i < left_up_location_ori.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(left_up_location_ori[i][0], left_up_location_ori[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	for (int i = 0; i < left_down_location_ori.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(left_down_location_ori[i][0], left_down_location_ori[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	for (int i = 0; i < right_up_location_ori.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(right_up_location_ori[i][0], right_up_location_ori[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	cv::imwrite("reOriImage.jpg", reOriImage);
}