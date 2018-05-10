#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <iostream>
#include <fstream> 
#include "fit.h"


using namespace std;


//����ͼ��Աȶ�
void ContrastImage(cv::Mat &image, int alpha, int beta);

//doubleת��Ϊstring
string doubleToString(double num);

//�ǵ�����ƶ�
vector<vector<double>> moveOri(vector<vector<double>> &o_location, vector<vector<double>> c_location, vector<vector<double>> eye_location);
//��Ļ�������������Ĵ洢���
void writeCsv(vector<vector<double>> ori_camera, string name);

//�������ͼ��
void writefor_p(cv::Mat image, vector<cv::Point> fitdata, int mode);


//��Ļͼ���Ԥ����
class location{

public:
	location();
	location(cv::Mat image);
	//ͨ������
	void splitImage();

	//��ȡ����
	void getCooder();

	//��ȡ���ĵ�
	void getCenter();

	//��ȡ�������
	void getFitdata();

	//��ȡ���
	void getMargin();

	//��������
	void sortCooder();

	//���޷���
	void separatefour();

	//������궨λ
	void fit_location();

	//ԭʼ���궨λ
	void ori_location();

	//�������ͼ��
	void writefitImage_help(cv::Mat image, vector<cv::Point> fitdata, int mode);

	//�������ͼ��
	void writeImage(cv::Mat image);

	//��дλ����Ϣ
	void write_message(cv::Mat image,int mode);

protected:
	int location::location_points(vector<vector<double>> &right_down_location, vector<cv::Point>right_down,
		cv::Point cameracenter_JZ, vector<double>RotaterowMargin_JZ, vector<double>RotatecolMargin_JZ, int Quadrant, int threshold);

private:

	//��ȡ�������
	vector<cv::Point> getFitData_help(vector<cv::Point> allCooder, cv::Point center_temp, cv::Point fitcenter_temp);

	//��ȡ����
	vector<cv::Point> getCooder_help(cv::Mat temp);

	//��������
	void sortCooder_help(vector<cv::Point>& pointes, int mode);

	//��ȡ������
	vector<double> getMargin_help(vector<cv::Point> &pointes, int mode);

	//��Ͻǵ㶨λ
	int find_mini_Point(vector<cv::Point> Pointes_temp, cv::Point center_temp);

	cv::Point locationAnddelete(vector<cv::Point> & right_down, int index, vector<vector<double>> &right_down_location,
		int row, int col, cv::Point center, vector<double> row_margin, vector<double> col_margin,
		int mode /*mode=1,ˮƽ�߲��ң�mode=2����ֱ�߲��ң�mode=3,б�߲���*/, int Quadrant/*����*/, int threshold);


	//ԭʼ���궨λ
	vector<vector<double>> ori_location_help(vector<cv::Point> Fit_othercooder, vector<vector<double>> fit_location, vector<cv::Point> othercooder, int&rownum, int &colnum);


	//��������
public:
	//ԭʼͼ�������
	vector<cv::Point> row_points;
	vector<cv::Point> col_points;
	vector<cv::Point> other_points;
	cv::Point center;

	//������
	vector<double> col_margin;
	vector<double> row_margin;

	//���޷���	
	vector<cv::Point> left_up, left_down, right_up, right_down;
	vector<vector<double>> left_up_location, left_down_location, right_up_location, right_down_location;
	vector<cv::Point> left_up_fit, left_down_fit, right_up_fit, right_down_fit;
	vector<vector<double>> left_up_fit_location, left_down_fit_location, right_up_fit_location, right_down_fit_location;

	//ͨ������
	cv::Mat ori_image;
	cv::Mat red_image;
	cv::Mat blue_image;
	cv::Mat green_image;

	//�������
	vector<cv::Point> row_points_fit;
	vector<cv::Point> col_points_fit;
	vector<cv::Point> other_points_fit;
	cv::Point center_fit;
	int rownum;
	int colnum;

};



//���ͼ��
class c_location :public location{

public:
	c_location() :location()
	{

	}
	c_location(cv::Mat image) :location(image)
	{
	}

	//ͨ������(�ض���)
	void splitImage();

	//����ϸ��
	void cooderRight();

	//���ͼ�����
	void undisPoints();

	//��ȡ������꣨�ض��壩
	void getFitdata();//���ֻ��Ҫ���и�ֵ��������

	//�ƶ���������ͼ��
	void movedata();

	//���ͼ����ת
	void CalrotatePoint();

	//������ͼ��λ ���ض��壩
	void fit_location();

	//���������Ĳο�ͼ��
	void generate_P_image();

private:
	//����ϸ��
	vector<cv::Point> cameraCooderRight(vector<cv::Point>camera_Redcooder, cv::Point cameraCenter, int mode);
	vector<cv::Point> otherCooderRight(const vector<cv::Point> camera_Redcooder_right,const vector<cv::Point> camera_Bluecooder_right, vector<cv::Point> camera_Greencooder);

	//���ͼ�����
	void UndistortPoints(cv::InputArray distorted, cv::OutputArray undistorted, cv::InputArray K, cv::InputArray D);
	void myUndistortPoints(std::vector<cv::Point> distorted, std::vector<cv::Point>& undistorted, cv::Mat K, cv::Mat D);

	//���ͼ����ת
	vector<cv::Point> CalrotatePoint_help(vector<cv::Point> points, cv::Point center, double angle);

	// �����������۾�ͼ��
	vector<vector<double>> c_location::generate_Eye_other(vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, int Quadrant, int rownum, int colnum);
	//����eyeͼ���row col
	void c_location::generate_row_col(int colnum, int rownum, vector<double> R_camerarow_JZ, vector<double> R_cameracol_JZ, cv::Point center, vector<cv::Point> &rowcooder_eye, vector<cv::Point>& colcooder_eye);


public:
	double angle;
	vector<vector<double>>  right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye;
};


cv::Mat reori_show(location temp);