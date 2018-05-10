#include"preprocessing.h"
#include"fit.h"


//增加图像对比度
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

string doubleToString(double num)
{
	return to_string(int(num));
}



//-------------------------- 屏幕图像 ------------------------//
location::location()
{
}
location::location(cv::Mat image)
{
	ori_image = image;
}

//分离三个通道
void location::splitImage()
{
	cv::Mat oriChannel[3];
	cv::split(ori_image, oriChannel);

	// red channel
	oriChannel[2].copyTo(red_image);

	// bule channel
	oriChannel[0].copyTo(blue_image);

	// green channel
	oriChannel[1].copyTo(green_image);
}

//获取坐标
vector<cv::Point> location::getCooder_help(cv::Mat temp)
{
	vector<cv::Point> colcooder;
	vector<vector<cv::Point>> contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(temp, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
	for (int i = 0; i < contours.size(); i++)
	{
		if (contours[i].size() < 5)
		{
			continue;
		}
		cv::RotatedRect box = cv::fitEllipse(contours[i]);
		colcooder.push_back(box.center);
	}
	return colcooder;
}
void location::getCooder()
{
	row_points = getCooder_help(blue_image);
	col_points = getCooder_help(red_image);
	other_points = getCooder_help(green_image);
}

//获取中心点坐标
void location::getCenter()
{
	vector<vector<cv::Point>> contours;
	vector<cv::Vec4i> hierarchy;

	cv::findContours(red_image, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0));
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
	if (contours[maxindex].size() >= 5)
	{
		cv::RotatedRect box = cv::fitEllipse(contours[maxindex]);
		center = box.center;
	}
}

//获取拟合数据
vector<cv::Point> location::getFitData_help(vector<cv::Point> allCooder, cv::Point center, cv::Point fitcenter)
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
void location::getFitdata()
{
	center_fit = cv::Point(900, 1000);
	//1. get fit_row cooder
	row_points_fit = getFitData_help(row_points, center, center_fit);

	//2. get fit_col cooder
	col_points_fit = getFitData_help(col_points, center, center_fit);

	//3. get fit_other cooder
	other_points_fit = getFitData_help(other_points, center, center_fit);
}

//坐标排序
void location::sortCooder_help(vector<cv::Point> &pointes, int mode)
{
	for (int i = 0; i < pointes.size(); i++)
	{
		for (int j = i; j < pointes.size(); j++)
		{
			if (mode == 0)
			{
				if (pointes[i].y > pointes[j].y)
				{
					cv::Point temp = pointes[i];
					pointes[i] = pointes[j];
					pointes[j] = temp;
				}
			}
			if (mode == 1)
			{
				if (pointes[i].x > pointes[j].x)
				{
					cv::Point temp = pointes[i];
					pointes[i] = pointes[j];
					pointes[j] = temp;
				}
			}
		}
	}
}
void location::sortCooder()
{
	sortCooder_help(row_points_fit, 1);
	sortCooder_help(col_points_fit, 0);
}

//获取间隔
vector<double> location::getMargin_help(vector<cv::Point>& pointes, int mode)
{
	vector<double> rowMargin;

	int centerNum = 0;
	for (int i = 0; i < pointes.size(); i++)
	{
		if (mode == 0)
		{
			if (abs(pointes[i].x - center_fit.x)<10)
			{
				centerNum = i;
			}
		}
		if (mode == 1)
		{
			if (abs(pointes[i].y - center_fit.y)<10)
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
			temp = pointes[ii].x - pointes[ii - 1].x;
		}
		if (mode == 1)
		{
			temp = temp = pointes[ii].y - pointes[ii - 1].y;
		}

		rowMargin.push_back(temp);
	}
	return rowMargin;
}
void location::getMargin()
{
	row_margin = getMargin_help(row_points_fit, 0);
	col_margin = getMargin_help(col_points_fit, 1);
}

//象限分离
void location::separatefour()
{
	for (int i = 0; i<other_points.size(); i++)
	{
		//left_up
		if (other_points[i].x<center.x && other_points[i].y<center.y)
		{
			left_up.push_back(other_points[i]);
		}
		//left_dowm
		if (other_points[i].x>center.x && other_points[i].y<center.y)
		{
			right_up.push_back(other_points[i]);
		}
		//right_up
		if (other_points[i].x<center.x && other_points[i].y>center.y)
		{
			left_down.push_back(other_points[i]);
		}
		//right_down
		if (other_points[i].x>center.x && other_points[i].y>center.y)
		{
			right_down.push_back(other_points[i]);
		}
	}


	for (int i = 0; i<other_points_fit.size(); i++)
	{
		//left_up
		if (other_points_fit[i].x<center_fit.x && other_points_fit[i].y<center_fit.y)
		{
			left_up_fit.push_back(other_points_fit[i]);
		}
		//left_dowm
		if (other_points_fit[i].x>center_fit.x && other_points_fit[i].y<center_fit.y)
		{
			right_up_fit.push_back(other_points_fit[i]);
		}
		//right_up
		if (other_points_fit[i].x<center_fit.x && other_points_fit[i].y>center_fit.y)
		{
			left_down_fit.push_back(other_points_fit[i]);
		}
		//right_down
		if (other_points_fit[i].x>center_fit.x && other_points_fit[i].y>center_fit.y)
		{
			right_down_fit.push_back(other_points_fit[i]);
		}
	}
}


//find mini Point
int location::find_mini_Point(vector<cv::Point> Pointes, cv::Point center)
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
cv::Point location::locationAnddelete(vector<cv::Point> & right_down_temp, int index, vector<vector<double>> &right_down_location, int row, int col, cv::Point center,
	vector<double> row_margin, vector<double> col_margin, int mode /*mode=1,水平边查找；mode=2，竖直边查找，mode=3,斜边查找*/, int Quadrant/*象限*/, int threshold)
{
	vector<double> location_temp;
	cv::Point rePoint;
	cv::Point miss_point;
	//判断是否有缺点
	double length = sqrt(pow(right_down_temp[index].x - center.x, 2) + pow(right_down_temp[index].y - center.y, 2));
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
			if (abs(right_down_temp[index].y - center.y) > 20)//水平方向上找点出现错误
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
				location_temp.push_back(right_down_temp[index].x);
				location_temp.push_back(right_down_temp[index].y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = right_down_temp[index];
				//删除找到的点
				right_down_temp.erase(right_down_temp.begin() + index);
			}

		}

		if (mode == 2)
		{
			if (abs(right_down_temp[index].x - center.x) > 20)//水平方向上找点出现错误
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
				location_temp.push_back(right_down_temp[index].x);
				location_temp.push_back(right_down_temp[index].y);
				location_temp.push_back(row);
				location_temp.push_back(col);
				rePoint = right_down_temp[index];
				//删除找到的点
				right_down_temp.erase(right_down_temp.begin() + index);
			}

		}

		if (mode == 3)
		{
			//right_down_location 接收
			location_temp.push_back(right_down_temp[index].x);
			location_temp.push_back(right_down_temp[index].y);
			location_temp.push_back(row);
			location_temp.push_back(col);
			rePoint = right_down_temp[index];
			//删除找到的点
			right_down_temp.erase(right_down_temp.begin() + index);
		}
	}
	right_down_location.push_back(location_temp);
	return rePoint;
}

//拟合坐标定位
int location::location_points(vector<vector<double>> &right_down_location, vector<cv::Point>right_down, cv::Point cameracenter_JZ, vector<double>RotaterowMargin_JZ, vector<double>RotatecolMargin_JZ,
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


void location::fit_location()
{
	location_points(right_up_fit_location, right_up_fit, center_fit, row_margin, col_margin, 1, 20);//第一象限
	location_points(left_up_fit_location, left_up_fit, center_fit, row_margin, col_margin, 2, 20);//第二象限
	location_points(left_down_fit_location, left_down_fit, center_fit, row_margin, col_margin, 3, 20);//第三象限
	location_points(right_down_fit_location, right_down_fit, center_fit, row_margin, col_margin, 4, 20);//第四象限
}

// 原始坐标定位
vector<vector<double>> location::ori_location_help(vector<cv::Point> Fit_othercooder, vector<vector<double>> fit_location, vector<cv::Point> othercooder, int&rownum, int &colnum)
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

//原始坐标定位
void location::ori_location()
{
	right_down_location = ori_location_help(other_points_fit, right_down_fit_location, other_points, rownum, colnum);
	left_up_location = ori_location_help(other_points_fit, left_up_fit_location, other_points, rownum, colnum);
	left_down_location = ori_location_help(other_points_fit, left_down_fit_location, other_points, rownum, colnum);
	right_up_location = ori_location_help(other_points_fit, right_up_fit_location, other_points, rownum, colnum);
}


//绘制拟合图像
void location::writefitImage_help(cv::Mat image, vector<cv::Point> fitdata, int mode)
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

//绘制拟合图像
void location::writeImage(cv::Mat image)
{
	writefitImage_help(image, row_points_fit, 1);
	writefitImage_help(image, col_points_fit, 0);
	writefitImage_help(image, other_points_fit, 2);
	cv::circle(image, center_fit, 15, cv::Scalar(255, 0, 255), -1, 8);

}

//编写位置信息
void location::write_message(cv::Mat image,int mode)
{
	int font_face = cv::FONT_HERSHEY_COMPLEX;

	if (mode == 1)
	{
		//写位置信息（坐标编码）
		for (int i = 0; i < right_down_fit_location.size(); i++)
		{
			string text = doubleToString(right_down_fit_location[i][2]).append(",").append(doubleToString(right_down_fit_location[i][3]));
			cv::putText(image, text, cv::Point(right_down_fit_location[i][0] - 5, right_down_fit_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < left_up_fit_location.size(); i++)
		{
			string text = doubleToString(left_up_fit_location[i][2]).append(",").append(doubleToString(left_up_fit_location[i][3]));
			cv::putText(image, text, cv::Point(left_up_fit_location[i][0] - 5, left_up_fit_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < left_down_fit_location.size(); i++)
		{
			string text = doubleToString(left_down_fit_location[i][2]).append(",").append(doubleToString(left_down_fit_location[i][3]));
			cv::putText(image, text, cv::Point(left_down_fit_location[i][0] - 5, left_down_fit_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < right_up_fit_location.size(); i++)
		{
			string text = doubleToString(right_up_fit_location[i][2]).append(",").append(doubleToString(right_up_fit_location[i][3]));
			cv::putText(image, text, cv::Point(right_up_fit_location[i][0] - 5, right_up_fit_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}

	}
	
	if (mode == 2)
	{
		//写位置信息（坐标编码）
		for (int i = 0; i < right_down_location.size(); i++)
		{
			string text = doubleToString(right_down_location[i][2]).append(",").append(doubleToString(right_down_location[i][3]));
			cv::putText(image, text, cv::Point(right_down_location[i][0] - 5, right_down_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < left_up_location.size(); i++)
		{
			string text = doubleToString(left_up_location[i][2]).append(",").append(doubleToString(left_up_location[i][3]));
			cv::putText(image, text, cv::Point(left_up_location[i][0] - 5, left_up_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < left_down_location.size(); i++)
		{
			string text = doubleToString(left_down_location[i][2]).append(",").append(doubleToString(left_down_location[i][3]));
			cv::putText(image, text, cv::Point(left_down_location[i][0] - 5, left_down_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
		for (int i = 0; i < right_up_location.size(); i++)
		{
			string text = doubleToString(right_up_location[i][2]).append(",").append(doubleToString(right_up_location[i][3]));
			cv::putText(image, text, cv::Point(right_up_location[i][0] - 5, right_up_location[i][1] + 8), font_face, 0.4, cv::Scalar(0, 0, 255), 1, 8);
		}
	}
}




//-------------------------- 相机图像 ------------------------//
//c_location 重定义splitImage函数
void c_location::splitImage()
{
	static cv::Mat element1 = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(15, 15));
	static cv::Mat element2 = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));

	cv::Mat cameraChannel[3];
	cv::split(ori_image, cameraChannel);
	//3. red channel
	cv::Mat cameraRed;
	cameraChannel[2].copyTo(cameraRed);

	cv::threshold(cameraRed, cameraRed, 100, 255, 0);
	//cv::imwrite("test.jpg", cameraRed);
	cv::erode(cameraRed, cameraRed, element1);
	cv::dilate(cameraRed, red_image, element2);

	//4. bule channel
	cv::Mat cameraBlue;
	cameraChannel[0].copyTo(cameraBlue);
	cv::threshold(cameraBlue, cameraBlue, 100, 255, 0);
	cv::erode(cameraBlue, cameraBlue, element1);
	cv::dilate(cameraBlue, blue_image, element2);

	//5. green channel 
	cv::Mat cameraGreen;
	cameraChannel[1].copyTo(cameraGreen);
	cv::threshold(cameraGreen, green_image, 100, 255, 0);

}

//c_location 特有的坐标细化
vector<cv::Point> c_location::cameraCooderRight(vector<cv::Point>camera_Redcooder, cv::Point cameraCenter, int mode)
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

vector<cv::Point> c_location::otherCooderRight(const vector<cv::Point> camera_Redcooder_right,
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

//坐标细化
void c_location::cooderRight()
{
	col_points = cameraCooderRight(col_points, center, 1);
	row_points = cameraCooderRight(row_points, center, 0);
	other_points = otherCooderRight(row_points, col_points, other_points);
}


//相机图像矫正
void c_location::UndistortPoints(cv::InputArray distorted, cv::OutputArray undistorted, cv::InputArray K, cv::InputArray D)
{
	cv::InputArray R = cv::noArray();
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

void c_location::myUndistortPoints(std::vector<cv::Point> distorted, std::vector<cv::Point>& undistorted, cv::Mat K, cv::Mat D)
{
	std::vector<cv::Point2f> tmpDistorted, tmpUndistorted;
	for (auto p : distorted)
		//for (int i = 0; i < distorted.size(); i++)
	{
		tmpDistorted.push_back(cv::Point2f((float)p.x, (float)p.y));
	}
	if (tmpDistorted.size() != distorted.size())
		std::cout << "point size error!" << std::endl;

	c_location::UndistortPoints(tmpDistorted, tmpUndistorted, K, D);

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

void c_location::undisPoints()
{
	cv::Mat intrinsicMatrix = cv::Mat::zeros(3, 3, CV_64F);
	cv::Mat distortionCoeff = cv::Mat::ones(1, 4, CV_64F);

	intrinsicMatrix.at<double>(0, 0) = 1601.5755521;
	intrinsicMatrix.at<double>(0, 1) = 0;
	intrinsicMatrix.at<double>(0, 2) = 2047.2178200;
	intrinsicMatrix.at<double>(1, 0) = 0;
	intrinsicMatrix.at<double>(1, 1) = 1602.774750;
	intrinsicMatrix.at<double>(1, 2) = 1481.77;
	intrinsicMatrix.at<double>(2, 0) = 0;
	intrinsicMatrix.at<double>(2, 1) = 0;
	intrinsicMatrix.at<double>(2, 2) = 1;

	distortionCoeff.at<double>(0, 0) = -0.00730083;
	distortionCoeff.at<double>(0, 1) = -0.000315583;
	distortionCoeff.at<double>(0, 2) = -0.000565603;
	distortionCoeff.at<double>(0, 3) = -0.000978041;


	/*cv::Mat intrinsicMatrix = cv::Mat::zeros(3, 3, CV_64F);
	cv::Mat distortionCoeff = cv::Mat::ones(1, 4, CV_64F);

	intrinsicMatrix.at<double>(0, 0) = 1622.74;
	intrinsicMatrix.at<double>(0, 1) = 0;
	intrinsicMatrix.at<double>(0, 2) = 2047.22;
	intrinsicMatrix.at<double>(1, 0) = 0;
	intrinsicMatrix.at<double>(1, 1) = 1624.77;
	intrinsicMatrix.at<double>(1, 2) = 1481.77;
	intrinsicMatrix.at<double>(2, 0) = 0;
	intrinsicMatrix.at<double>(2, 1) = 0;
	intrinsicMatrix.at<double>(2, 2) = 1;

	distortionCoeff.at<double>(0, 0) = -0.0103251;
	distortionCoeff.at<double>(0, 1) = -0.00183935;
	distortionCoeff.at<double>(0, 2) = -0.00103306;
	distortionCoeff.at<double>(0, 3) = -0.00116156;*/


	vector<cv::Point> cameracenter_right;
	cameracenter_right.push_back(center);

	vector<cv::Point> cameracenter_JZ_temp;
	myUndistortPoints(cameracenter_right, cameracenter_JZ_temp, intrinsicMatrix, distortionCoeff);
	center = cameracenter_JZ_temp[0];

	//get JZ image row cooder
	vector<cv::Point> row_temp;
	myUndistortPoints(row_points, row_temp, intrinsicMatrix, distortionCoeff);
	row_points = row_temp;

	//get JZ image col cooder
	vector<cv::Point> col_temp;
	myUndistortPoints(col_points, col_temp, intrinsicMatrix, distortionCoeff);
	col_points = col_temp;

	// get JZ image other cooder
	vector<cv::Point> other_temp;
	myUndistortPoints(other_points, other_temp, intrinsicMatrix, distortionCoeff);
	other_points = other_temp;
}


//获取拟合坐标
void c_location::getFitdata()
{
	row_points_fit = row_points;
	col_points_fit = col_points;
	other_points_fit = other_points;
	center_fit = center;
}

void c_location::movedata()
{
	//jcameracenter_JZ camerarow_JZ,cameracol_JZ,cameraother_JZ 全部平移
	for (int i = 0; i < row_points.size(); i++)
	{
		row_points[i].x += 2500 - center.x;
		row_points[i].y += 2000 - center.y;
	}
	for (int i = 0; i < col_points.size(); i++)
	{
		col_points[i].x += 2500 - center.x;
		col_points[i].y += 2000 - center.y;
	}
	for (int i = 0; i < other_points.size(); i++)
	{
		other_points[i].x += 2500 - center.x;
		other_points[i].y += 2000 - center.y;
	}
	center.x += 2500 - center.x;
	center.y += 2000 - center.y;
}

//相机图像旋转
// calclate the cooder of rotate
vector<cv::Point> c_location::CalrotatePoint_help(vector<cv::Point> points, cv::Point center, double angle)
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

void c_location::CalrotatePoint()
{
	// rotate JZ image
	double angle1, angle2, angle3, angle4;

	double temp_left = atan2((row_points_fit[0].y - center_fit.y), (center_fit.x - row_points_fit[0].x));
	angle1 = temp_left;

	double temp_up = atan2((col_points_fit[0].x - center_fit.x), (center_fit.y - col_points_fit[0].y));
	angle2 = -temp_up;

	double temp_right = atan2((row_points_fit[row_points_fit.size() - 1].y - center_fit.y), (row_points_fit[row_points_fit.size() - 1].x - center_fit.x));
	angle3 = -temp_right;

	double temp_down = atan2((col_points_fit[col_points_fit.size() - 1].x - center_fit.x), (col_points_fit[col_points_fit.size() - 1].y - center_fit.y));
	angle4 = temp_down;

	angle = -(angle1 + angle2 + angle3 + angle4) / 4;


	row_points_fit = CalrotatePoint_help(row_points_fit, center_fit, angle);

	col_points_fit = CalrotatePoint_help(col_points_fit, center_fit, angle);

	other_points_fit = CalrotatePoint_help(other_points_fit, center_fit, angle);

}


//相机拟合图像定位（重定义）
void c_location::fit_location()
{
	vector<cv::Point> right_up_temp = right_up;
	location_points(right_up_fit_location, right_up_fit, center_fit, row_margin, col_margin, 1, 61);//第一象限
	location_points(left_up_fit_location, left_up_fit, center_fit, row_margin, col_margin, 2, 61);//第二象限
	location_points(left_down_fit_location, left_down_fit, center_fit, row_margin, col_margin, 3, 61);//第三象限
	location_points(right_down_fit_location, right_down_fit, center_fit, row_margin, col_margin, 4, 61);//第四象限
}


//生成完美图像
// 生成完美的眼睛图像
vector<vector<double>> c_location::generate_Eye_other(vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, int Quadrant, int rownum, int colnum)
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
void c_location::generate_row_col(int colnum, int rownum, vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, vector<cv::Point> &rowcooder_eye, vector<cv::Point>& colcooder_eye)
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
		if (y_index < 0)
		{
			temp1.y = home1.y - R_cameracol_JZ[y_index] - 10;
		}
		else
		{
			temp1.y = home1.y - R_cameracol_JZ[y_index];
		}
		home1.y = temp1.y;

		temp2.x = center.x;
		if (y_index < 0)
		{
			temp2.y = home2.y + R_cameracol_JZ[y_index] - 10;
		}
		else
		{
			temp2.y = home2.y + R_cameracol_JZ[y_index];
		}
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
		if (x_index < 0)
		{
			temp1.x = home1.x - R_camerarow_JZ[x_index];
		}
		else
		{
			temp1.x = home1.x - R_camerarow_JZ[x_index];
		}
		home1.x = temp1.x;

		cv::Point temp2;
		temp2.y = center.y;
		if (x_index < 0)
		{
			temp2.x = home2.x + R_camerarow_JZ[x_index];
		}
		else
		{
			temp2.x = home2.x + R_camerarow_JZ[x_index];
		}
		//temp2.x = home2.x + R_camerarow_JZ[x_index];
		home2.x = temp2.x;

		rowcooder_eye.push_back(temp1);
		rowcooder_eye.push_back(temp2);
	}
	colcooder_eye.push_back(center);
	rowcooder_eye.push_back(center);

}


void c_location::generate_P_image()
{
	cv::Mat prefectEyeImage(4000, 5000, CV_8UC3, cv::Scalar::all(0));

	//生成prefecteye 图像的参数；
	vector<cv::Point> rowcooder_eye, colcooder_eye;
	generate_row_col(colnum, rownum, row_margin, col_margin, center, rowcooder_eye, colcooder_eye);
	cv::Point center_eye = center;
	writefitImage_help(prefectEyeImage, rowcooder_eye, 1);
	writefitImage_help(prefectEyeImage, colcooder_eye, 0);
	cv::circle(prefectEyeImage, center_eye, 15, cv::Scalar(255, 0, 255), -1, 8);

	
	//四部分点的位置信息
	right_up_location_Eye = generate_Eye_other(row_margin, col_margin, center, 1, rownum, colnum);
	left_up_location_Eye = generate_Eye_other(row_margin, col_margin, center, 2, rownum, colnum);
	left_down_location_Eye = generate_Eye_other(row_margin, col_margin, center, 3, rownum, colnum);
	right_down_location_Eye = generate_Eye_other(row_margin, col_margin, center, 4, rownum, colnum);

	int font_face = cv::FONT_HERSHEY_COMPLEX;
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




//角点迭代移动
vector<vector<double>> moveOri(vector<vector<double>> &o_location, vector<vector<double>> c_location, vector<vector<double>> eye_location)
{
	//记录参与调整的点的个数
	int num_adjust = 0;
	int num_adjust_no = 0;

	////存储参与调整之后的屏幕坐标及相机图像的坐标
	vector<vector<double>> ori_camera_location;

	//存储在c_location 中的位置索引，及eye_location中的位置索引
	int cam_index = -1;
	int eye_index = -1;
	for (int i = 0; i < o_location.size(); i++)
	{


		//查找相机位置信息中与屏幕图像位置信息一致的点
		for (int j = 0; j < c_location.size(); j++)
		{
			if (o_location[i][2] == c_location[j][2] &&
				o_location[i][3] == c_location[j][3])
			{
				cam_index = j;
				break;
			}
		}

		//查找eye图像中与屏幕位置信息一致的点
		for (int k = 0; k < eye_location.size(); k++)
		{
			if (o_location[i][2] == eye_location[k][2] &&
				o_location[i][3] == eye_location[k][3])
			{
				eye_index = k;
				break;
			}
		}

		if (eye_index == -1 || cam_index == -1)
		{
			cout << "角点对应时：：：没找到对应的点" << endl;
		}
		else
		{
			
			//找的对应点进行移动修正
			//判断camImage上的点与eyeimage上的点到底相差多少
			double x_length = c_location[cam_index][0] - eye_location[eye_index][0];
			double y_length = c_location[cam_index][1] - eye_location[eye_index][1];

			if (sqrt(pow(x_length, 2) + pow(y_length, 2))<4)
			{
				num_adjust_no += 1;
				
			}
			else
			{
				num_adjust += 1;

				if (x_length<0)
				{
					o_location[i][0] += 1;
				}
				else
				{
					o_location[i][0] -= 1;
				}

				if (y_length<0)
				{
					o_location[i][1] += 1;
				}
				else
				{
					o_location[i][1] -= 1;
				}
			}

			vector<double> temp;
			temp.push_back(o_location[i][0]);
			temp.push_back(o_location[i][1]);
			temp.push_back(c_location[cam_index][0]);
			temp.push_back(c_location[cam_index][1]);
			temp.push_back(eye_location[eye_index][0]);
			temp.push_back(eye_location[eye_index][1]);
			//保存ori坐标与camera坐标
			ori_camera_location.push_back(temp);
		}
	}

	cout << "参与调整 " << to_string(num_adjust) << endl;
	cout << "没有参与调整 " << to_string(num_adjust_no) << endl;

	return ori_camera_location;
}

//屏幕坐标与相机坐标的存储输出
void writeCsv(vector<vector<double>> ori_camera, string name)
{
	//保存ori坐标与camera坐标到csv文件
	ofstream ofile;
	ofile.open(name, ios::out | ios::trunc); //判断.csv文件是否存在，不存在则创建
	ofile << "ori-x" << "," << "ori-y" << "," << "camera-x" << "," << "camera-y" << "," << "prefrct-x" << "," << "prefrct-y" << endl;

	for (int i = 0; i < ori_camera.size(); i++)
	{
		ofile << ori_camera[i][0] << "," << ori_camera[i][1] << ","
			<< ori_camera[i][2] << "," << ori_camera[i][3] << "," << ori_camera[i][4] << "," << ori_camera[i][5] << endl;
	}
	ofile.close();
}

//绘制调整后的屏幕图像
void writefor_p(cv::Mat image, vector<cv::Point> fitdata, int mode)
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

cv::Mat reori_show(location temp)
{
	//重新生成原始图片
	cv::Mat reOriImage(temp.ori_image.size(), CV_8UC3, cv::Scalar::all(0));
	writefor_p(reOriImage, temp.row_points, 1);
	writefor_p(reOriImage, temp.col_points, 0);
	cv::circle(reOriImage, temp.center, 15, cv::Scalar(255, 0, 255), -1, 8);

	for (int i = 0; i < temp.right_down_location.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(temp.right_down_location[i][0], temp.right_down_location[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	for (int i = 0; i < temp.left_up_location.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(temp.left_up_location[i][0], temp.left_up_location[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	for (int i = 0; i < temp.left_down_location.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(temp.left_down_location[i][0], temp.left_down_location[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	for (int i = 0; i < temp.right_up_location.size(); i++)
	{
		cv::circle(reOriImage, cv::Point(temp.right_up_location[i][0], temp.right_up_location[i][1]), 7, cv::Scalar(0, 255, 0), -1, 8);
	}
	cv::imwrite("reOriImage.jpg", reOriImage);

	return reOriImage;
}