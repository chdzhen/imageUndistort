#include "preprocessing.h"

int ttt_main()
{

	//---------------------屏幕图像预处理---------------------//
	//读取屏幕图像
	cv::Mat image;
	image = cv::imread("ori.bmp");
	if (image.data == NULL)
	{
		cout << "ori_image is null" << endl;
		return 1;
	}

	//初始化屏幕图像对象
	location ori_location(image);

	//通道分离
	ori_location.splitImage();

	//获取角点坐标
	ori_location.getCooder();
	ori_location.getCenter();

	//获取拟合图像
	ori_location.getFitdata();

	//------测试：：拟合图像绘制
	cv::Mat fitImage(1800, 2000, CV_8UC3, cv::Scalar::all(0));
	ori_location.writeImage(fitImage);
	cv::imwrite("fitimage.jpg", fitImage);

	//坐标排序
	ori_location.sortCooder();

	//获取间隔
	ori_location.getMargin();

	//象限分离
	ori_location.separatefour();

	//拟合坐标定位
	ori_location.fit_location();

	//------测试：：拟合图像位置信息编写
	ori_location.write_message(fitImage, 1);
	cv::imwrite("fitimage.jpg", fitImage);

	//屏幕图像定位
	ori_location.ori_location();

	//------测试：：屏幕图像位置信息编写
	cv::Mat oriLocation = cv::imread("oriLocation.bmp");
	ori_location.write_message(oriLocation, 2);
	cv::imwrite("oriLocation.bmp", oriLocation);


	for (int num = 0;; num++)
	{
		//---------------------相机图像预处理---------------------//

		//读取相机图像
		cv::Mat cameraImage;
		cameraImage = cv::imread("cameraImage.bmp");
		if (cameraImage.data == NULL)
		{
			cout << "cameraImage is null" << endl;
			return 2;
		}

		//初始化相机图像对象
		c_location camera_Location(cameraImage);

		//通道分离
		camera_Location.splitImage();

		//------测试：：显示分离的通道
		//cv::imwrite("test.bmp", camera_Location.green_image);
		//cv::waitKey();

		//获取角点坐标----重定义
		camera_Location.getCooder();
		camera_Location.getCenter();

		//坐标细化----特有的
		camera_Location.cooderRight();

		//相机图像矫正
		camera_Location.undisPoints();

		//移动拟合图像
		camera_Location.movedata();

		//移动标定的中心点
		cv::Point cli_center = cv::Point(2047.22, 1481.77);
		cli_center.x += 2500 - camera_Location.center.x;
		cli_center.y += 2000 - camera_Location.center.y;

		//获取拟合图像
		camera_Location.getFitdata();

		//相机图像矫正完之后，是相机图像的原始图像了，通过相机图像的拟合图像进行定位，此时拟合坐标已经进行了赋值
		cv::Mat cameraImage_location(4000, 5000, CV_8SC3, cv::Scalar::all(0));
		camera_Location.writeImage(cameraImage_location);
		cv::imwrite("cameraImage_location.bmp", cameraImage_location);



		//------测试：：拟合图像绘制
		cv::Mat cameraImage_fit(4000, 5000, CV_8UC3, cv::Scalar::all(0));
		camera_Location.writeImage(cameraImage_fit);
		cv::imwrite("cameraImage_fit.bmp", cameraImage_fit);

		//坐标排序
		camera_Location.sortCooder();

		//相机图像旋转
		camera_Location.CalrotatePoint();

		//获取间隔
		camera_Location.getMargin();

		//象限分离
		camera_Location.separatefour();

		//拟合坐标定位----(重定义)
		camera_Location.fit_location();

		//------测试：：拟合图像位置信息编写
		camera_Location.write_message(cameraImage_fit, 1);
		cv::imwrite("cameraImage_fit.bmp", cameraImage_fit);

		//屏幕图像定位
		camera_Location.ori_location();

		//------测试：：相机图像位置信息编写
		camera_Location.write_message(cameraImage_location, 2);
		cv::imwrite("cameraImage_location.bmp", cameraImage_location);

		//生成完美的图像
		camera_Location.generate_P_image();


		//----------移动迭代模块
		//根据camimage图像的坐标调整eyeimage的坐标

		//保存ori坐标与camera坐标
		vector<vector<double>> right_down_ori_camera;
		vector<vector<double>> left_up_ori_camera;
		vector<vector<double>> left_down_ori_camera;
		vector<vector<double>> right_up_ori_camera;

		//调整原始图片
		right_down_ori_camera = moveOri(ori_location.right_down_location, camera_Location.right_down_location, camera_Location.right_down_location_Eye);
		left_up_ori_camera = moveOri(ori_location.left_up_location, camera_Location.left_up_location, camera_Location.left_up_location_Eye);
		left_down_ori_camera = moveOri(ori_location.left_down_location, camera_Location.left_down_location, camera_Location.left_down_location_Eye);
		right_up_ori_camera = moveOri(ori_location.right_up_location, camera_Location.right_up_location, camera_Location.right_up_location_Eye);

		//重新生成屏幕图像
		cv::Mat reOriImage = reori_show(ori_location);

		//存储对应的角点坐标
		writeCsv(right_down_ori_camera, "right_down_ori_camera.csv");
		writeCsv(left_up_ori_camera, "left_up_ori_camera.csv");
		writeCsv(left_down_ori_camera, "left_down_ori_camera.csv");
		writeCsv(right_up_ori_camera, "right_up_ori_camera.csv");

		cout << "标记的中心点坐标为：";
		cout << cli_center.x << endl;
		cout << cli_center.y << endl;

		cout << "检测的中心点的坐标为：";
		cout << camera_Location.center.x << endl;
		cout << camera_Location.center.y << endl;

	}







	return 0;
}