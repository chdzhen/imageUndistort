#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <iostream>
#include <fstream> 
#include "fit.h"


using namespace std;


//增加图像对比度
void ContrastImage(cv::Mat &image, int alpha, int beta);

//double转换为string
string doubleToString(double num);

//角点迭代移动
vector<vector<double>> moveOri(vector<vector<double>> &o_location, vector<vector<double>> c_location, vector<vector<double>> eye_location);
//屏幕坐标与相机坐标的存储输出
void writeCsv(vector<vector<double>> ori_camera, string name);

//绘制拟合图像
void writefor_p(cv::Mat image, vector<cv::Point> fitdata, int mode);


//屏幕图像的预处理
class location{

public:
	location();
	location(cv::Mat image);
	//通道分离
	void splitImage();

	//获取坐标
	void getCooder();

	//获取中心点
	void getCenter();

	//获取拟合坐标
	void getFitdata();

	//获取间隔
	void getMargin();

	//坐标排序
	void sortCooder();

	//象限分离
	void separatefour();

	//拟合坐标定位
	void fit_location();

	//原始坐标定位
	void ori_location();

	//绘制拟合图像
	void writefitImage_help(cv::Mat image, vector<cv::Point> fitdata, int mode);

	//绘制拟合图像
	void writeImage(cv::Mat image);

	//编写位置信息
	void write_message(cv::Mat image,int mode);

protected:
	int location::location_points(vector<vector<double>> &right_down_location, vector<cv::Point>right_down,
		cv::Point cameracenter_JZ, vector<double>RotaterowMargin_JZ, vector<double>RotatecolMargin_JZ, int Quadrant, int threshold);

private:

	//获取拟合坐标
	vector<cv::Point> getFitData_help(vector<cv::Point> allCooder, cv::Point center_temp, cv::Point fitcenter_temp);

	//获取坐标
	vector<cv::Point> getCooder_help(cv::Mat temp);

	//坐标排序
	void sortCooder_help(vector<cv::Point>& pointes, int mode);

	//获取坐标间隔
	vector<double> getMargin_help(vector<cv::Point> &pointes, int mode);

	//拟合角点定位
	int find_mini_Point(vector<cv::Point> Pointes_temp, cv::Point center_temp);

	cv::Point locationAnddelete(vector<cv::Point> & right_down, int index, vector<vector<double>> &right_down_location,
		int row, int col, cv::Point center, vector<double> row_margin, vector<double> col_margin,
		int mode /*mode=1,水平边查找；mode=2，竖直边查找，mode=3,斜边查找*/, int Quadrant/*象限*/, int threshold);


	//原始坐标定位
	vector<vector<double>> ori_location_help(vector<cv::Point> Fit_othercooder, vector<vector<double>> fit_location, vector<cv::Point> othercooder, int&rownum, int &colnum);


	//公共属性
public:
	//原始图像的坐标
	vector<cv::Point> row_points;
	vector<cv::Point> col_points;
	vector<cv::Point> other_points;
	cv::Point center;

	//坐标间隔
	vector<double> col_margin;
	vector<double> row_margin;

	//象限分离	
	vector<cv::Point> left_up, left_down, right_up, right_down;
	vector<vector<double>> left_up_location, left_down_location, right_up_location, right_down_location;
	vector<cv::Point> left_up_fit, left_down_fit, right_up_fit, right_down_fit;
	vector<vector<double>> left_up_fit_location, left_down_fit_location, right_up_fit_location, right_down_fit_location;

	//通道分离
	cv::Mat ori_image;
	cv::Mat red_image;
	cv::Mat blue_image;
	cv::Mat green_image;

	//拟合坐标
	vector<cv::Point> row_points_fit;
	vector<cv::Point> col_points_fit;
	vector<cv::Point> other_points_fit;
	cv::Point center_fit;
	int rownum;
	int colnum;

};



//相机图像
class c_location :public location{

public:
	c_location() :location()
	{

	}
	c_location(cv::Mat image) :location(image)
	{
	}

	//通道分离(重定义)
	void splitImage();

	//坐标细化
	void cooderRight();

	//相机图像矫正
	void undisPoints();

	//获取拟合坐标（重定义）
	void getFitdata();//这个只需要进行赋值操作即可

	//移动相机的拟合图像
	void movedata();

	//相机图像旋转
	void CalrotatePoint();

	//相机拟合图像定位 （重定义）
	void fit_location();

	//生成完美的参考图像
	void generate_P_image();

private:
	//坐标细化
	vector<cv::Point> cameraCooderRight(vector<cv::Point>camera_Redcooder, cv::Point cameraCenter, int mode);
	vector<cv::Point> otherCooderRight(const vector<cv::Point> camera_Redcooder_right,const vector<cv::Point> camera_Bluecooder_right, vector<cv::Point> camera_Greencooder);

	//相机图像矫正
	void UndistortPoints(cv::InputArray distorted, cv::OutputArray undistorted, cv::InputArray K, cv::InputArray D);
	void myUndistortPoints(std::vector<cv::Point> distorted, std::vector<cv::Point>& undistorted, cv::Mat K, cv::Mat D);

	//相机图像旋转
	vector<cv::Point> CalrotatePoint_help(vector<cv::Point> points, cv::Point center, double angle);

	// 生成完美的眼睛图像
	vector<vector<double>> c_location::generate_Eye_other(vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, int Quadrant, int rownum, int colnum);
	//生成eye图像的row col
	void c_location::generate_row_col(int colnum, int rownum, vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, vector<cv::Point> &rowcooder_eye, vector<cv::Point>& colcooder_eye);


public:
	double angle;
	vector<vector<double>>  right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye;
};


cv::Mat reori_show(location temp);