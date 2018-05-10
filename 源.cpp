#include "preprocessing.h"

int ttt_main()
{

	//---------------------��Ļͼ��Ԥ����---------------------//
	//��ȡ��Ļͼ��
	cv::Mat image;
	image = cv::imread("ori.bmp");
	if (image.data == NULL)
	{
		cout << "ori_image is null" << endl;
		return 1;
	}

	//��ʼ����Ļͼ�����
	location ori_location(image);

	//ͨ������
	ori_location.splitImage();

	//��ȡ�ǵ�����
	ori_location.getCooder();
	ori_location.getCenter();

	//��ȡ���ͼ��
	ori_location.getFitdata();

	//------���ԣ������ͼ�����
	cv::Mat fitImage(1800, 2000, CV_8UC3, cv::Scalar::all(0));
	ori_location.writeImage(fitImage);
	cv::imwrite("fitimage.jpg", fitImage);

	//��������
	ori_location.sortCooder();

	//��ȡ���
	ori_location.getMargin();

	//���޷���
	ori_location.separatefour();

	//������궨λ
	ori_location.fit_location();

	//------���ԣ������ͼ��λ����Ϣ��д
	ori_location.write_message(fitImage, 1);
	cv::imwrite("fitimage.jpg", fitImage);

	//��Ļͼ��λ
	ori_location.ori_location();

	//------���ԣ�����Ļͼ��λ����Ϣ��д
	cv::Mat oriLocation = cv::imread("oriLocation.bmp");
	ori_location.write_message(oriLocation, 2);
	cv::imwrite("oriLocation.bmp", oriLocation);


	for (int num = 0;; num++)
	{
		//---------------------���ͼ��Ԥ����---------------------//

		//��ȡ���ͼ��
		cv::Mat cameraImage;
		cameraImage = cv::imread("cameraImage.bmp");
		if (cameraImage.data == NULL)
		{
			cout << "cameraImage is null" << endl;
			return 2;
		}

		//��ʼ�����ͼ�����
		c_location camera_Location(cameraImage);

		//ͨ������
		camera_Location.splitImage();

		//------���ԣ�����ʾ�����ͨ��
		//cv::imwrite("test.bmp", camera_Location.green_image);
		//cv::waitKey();

		//��ȡ�ǵ�����----�ض���
		camera_Location.getCooder();
		camera_Location.getCenter();

		//����ϸ��----���е�
		camera_Location.cooderRight();

		//���ͼ�����
		camera_Location.undisPoints();

		//�ƶ����ͼ��
		camera_Location.movedata();

		//�ƶ��궨�����ĵ�
		cv::Point cli_center = cv::Point(2047.22, 1481.77);
		cli_center.x += 2500 - camera_Location.center.x;
		cli_center.y += 2000 - camera_Location.center.y;

		//��ȡ���ͼ��
		camera_Location.getFitdata();

		//���ͼ�������֮�������ͼ���ԭʼͼ���ˣ�ͨ�����ͼ������ͼ����ж�λ����ʱ��������Ѿ������˸�ֵ
		cv::Mat cameraImage_location(4000, 5000, CV_8SC3, cv::Scalar::all(0));
		camera_Location.writeImage(cameraImage_location);
		cv::imwrite("cameraImage_location.bmp", cameraImage_location);



		//------���ԣ������ͼ�����
		cv::Mat cameraImage_fit(4000, 5000, CV_8UC3, cv::Scalar::all(0));
		camera_Location.writeImage(cameraImage_fit);
		cv::imwrite("cameraImage_fit.bmp", cameraImage_fit);

		//��������
		camera_Location.sortCooder();

		//���ͼ����ת
		camera_Location.CalrotatePoint();

		//��ȡ���
		camera_Location.getMargin();

		//���޷���
		camera_Location.separatefour();

		//������궨λ----(�ض���)
		camera_Location.fit_location();

		//------���ԣ������ͼ��λ����Ϣ��д
		camera_Location.write_message(cameraImage_fit, 1);
		cv::imwrite("cameraImage_fit.bmp", cameraImage_fit);

		//��Ļͼ��λ
		camera_Location.ori_location();

		//------���ԣ������ͼ��λ����Ϣ��д
		camera_Location.write_message(cameraImage_location, 2);
		cv::imwrite("cameraImage_location.bmp", cameraImage_location);

		//����������ͼ��
		camera_Location.generate_P_image();


		//----------�ƶ�����ģ��
		//����camimageͼ����������eyeimage������

		//����ori������camera����
		vector<vector<double>> right_down_ori_camera;
		vector<vector<double>> left_up_ori_camera;
		vector<vector<double>> left_down_ori_camera;
		vector<vector<double>> right_up_ori_camera;

		//����ԭʼͼƬ
		right_down_ori_camera = moveOri(ori_location.right_down_location, camera_Location.right_down_location, camera_Location.right_down_location_Eye);
		left_up_ori_camera = moveOri(ori_location.left_up_location, camera_Location.left_up_location, camera_Location.left_up_location_Eye);
		left_down_ori_camera = moveOri(ori_location.left_down_location, camera_Location.left_down_location, camera_Location.left_down_location_Eye);
		right_up_ori_camera = moveOri(ori_location.right_up_location, camera_Location.right_up_location, camera_Location.right_up_location_Eye);

		//����������Ļͼ��
		cv::Mat reOriImage = reori_show(ori_location);

		//�洢��Ӧ�Ľǵ�����
		writeCsv(right_down_ori_camera, "right_down_ori_camera.csv");
		writeCsv(left_up_ori_camera, "left_up_ori_camera.csv");
		writeCsv(left_down_ori_camera, "left_down_ori_camera.csv");
		writeCsv(right_up_ori_camera, "right_up_ori_camera.csv");

		cout << "��ǵ����ĵ�����Ϊ��";
		cout << cli_center.x << endl;
		cout << cli_center.y << endl;

		cout << "�������ĵ������Ϊ��";
		cout << camera_Location.center.x << endl;
		cout << camera_Location.center.y << endl;

	}







	return 0;
}