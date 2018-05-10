
#include <PvSampleUtils.h>
#include <PvBuffer.h>
#include <PvPixelType.h>
#include <PvBufferWriter.h>
#include <PvString.h>
#include <string>
#include <vector>
#include<io.h>
#include <fstream>
#include "preprocessing.h"



using namespace std;

PV_INIT_SIGNAL_HANDLER();




#include <PvSampleUtils.h>
#include <PvDevice.h>
#include <PvBuffer.h>
#include <PvStream.h>
#ifdef PV_GUI_NOT_AVAILABLE
#include <PvSystem.h>
#else
#include <PvDeviceFinderWnd.h>
#endif // PV_GUI_NOT_AVAILABLE

PV_INIT_SIGNAL_HANDLER();

#define BUFFER_COUNT ( 16 )

cv::Mat cameraImage = cv::Mat(3000, 4000, CV_8UC3, cv::Scalar(0, 0, 0));;


bool AcquireImages(cv::Mat oriImage)
{
	PvDeviceInfo *lDeviceInfo = NULL;
	PvResult lResult;

#ifdef PV_GUI_NOT_AVAILABLE
	PvSystem lSystem;
	lDeviceInfo = PvSelectDevice(lSystem);
	if (lDeviceInfo == NULL)
	{
		cout << "No device selected." << endl;
		return false;
	}
#else
	PvDeviceFinderWnd lDeviceFinderWnd;
	lResult = lDeviceFinderWnd.ShowModal();

	if (!lResult.IsOK())
	{
		cout << "No device selected." << endl;
		return false;
	}

	lDeviceInfo = lDeviceFinderWnd.GetSelected();
#endif // PV_GUI_NOT_AVAILABLE

	// Connect to the GEV Device
	PvDevice lDevice;
	cout << "Connecting to " << lDeviceInfo->GetMACAddress().GetAscii() << endl;

	if (!lDevice.Connect(lDeviceInfo).IsOK())
	{
		cout << "Unable to connect to " << lDeviceInfo->GetMACAddress().GetAscii() << endl;
		return false;
	}
	cout << "Successfully connected to " << lDeviceInfo->GetMACAddress().GetAscii() << endl << endl;

	// Get device parameters need to control streaming
	PvGenParameterArray *lDeviceParams = lDevice.GetGenParameters();
	PvGenInteger *lTLLocked = dynamic_cast<PvGenInteger *>(lDeviceParams->Get("TLParamsLocked"));
	PvGenInteger *lPayloadSize = dynamic_cast<PvGenInteger *>(lDeviceParams->Get("PayloadSize"));
	PvGenCommand *lStart = dynamic_cast<PvGenCommand *>(lDeviceParams->Get("AcquisitionStart"));
	PvGenCommand *lStop = dynamic_cast<PvGenCommand *>(lDeviceParams->Get("AcquisitionStop"));

	// Negotiate streaming packet size
	lDevice.NegotiatePacketSize();

	// Create the PvStream object
	PvStream lStream;

	// Open stream - have the PvDevice do it for us
	cout << "Opening stream to device" << endl;
	lStream.Open(lDeviceInfo->GetIPAddress());

	// Reading payload size from device
	PvInt64 lSize = 0;
	lPayloadSize->GetValue(lSize);

	// Use min of BUFFER_COUNT and how many buffers can be queued in PvStream
	PvUInt32 lBufferCount = (lStream.GetQueuedBufferMaximum() < BUFFER_COUNT) ?
		lStream.GetQueuedBufferMaximum() :
		BUFFER_COUNT;

	// Create, alloc buffers
	PvBuffer *lBuffers = new PvBuffer[lBufferCount];
	for (PvUInt32 i = 0; i < lBufferCount; i++)
	{
		lBuffers[i].Alloc(static_cast<PvUInt32>(lSize));
	}

	// Have to set the Device IP destination to the Stream
	lDevice.SetStreamDestination(lStream.GetLocalIPAddress(), lStream.GetLocalPort());


	// Queue all buffers in the stream
	for (PvUInt32 i = 0; i < lBufferCount; i++)
	{
		lStream.QueueBuffer(lBuffers + i);
	}

	// TLParamsLocked is optional but when present, it MUST be set to 1
	// before sending the AcquisitionStart command
	if (lTLLocked != NULL)
	{
		cout << "Setting TLParamsLocked to 1" << endl;
		lTLLocked->SetValue(1);
	}

	cout << "Resetting timestamp counter..." << endl;
	PvGenCommand *lResetTimestamp = dynamic_cast<PvGenCommand *>(lDeviceParams->Get("GevTimestampControlReset"));
	lResetTimestamp->Execute();

	// The buffers are queued in the stream, we just have to tell the device
	// to start sending us images
	cout << "Sending StartAcquisition command to device" << endl;
	lResult = lStart->Execute();



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

	//存储完美图像的坐标及位置信息
	vector<vector<double>> right_up, right_down, left_up, left_down;


	int num_loop = 0;

	int  index_camera = 0;

	while (!PvKbHit())
	{
		index_camera += 1;
		PvBuffer *lPtr = NULL;
		PvResult lOperationResult;

		// Retrieve next buffer		
		PvResult lResult = lStream.RetrieveBuffer(&lPtr, &lOperationResult, 1000);
		if (lResult.IsOK())
		{
			if (lOperationResult.IsOK())
			{
				// If the buffer contains an image, display width and height
				PvUInt32 lWidth = 0, lHeight = 0;
				if (lPtr->GetPayloadType() == PvPayloadTypeImage)
				{
					// Read width, height
					lWidth = lPtr->GetImage()->GetWidth();
					lHeight = lPtr->GetImage()->GetHeight();

					PvBuffer * lBufferRGB32 = new PvBuffer();
					lBufferRGB32->GetImage()->Alloc(lWidth, lHeight, PvPixelBGR8);

					PvBufferConverter lBufferConverter;
					lResult = lBufferConverter.Convert(lPtr, lBufferRGB32, true);

					cout << "num: " << index_camera << endl;

					PvUInt8 *lBuffrtPointr = lBufferRGB32->GetDataPointer();

					if (index_camera % 20 == 0)
					{
						for (int i = 0; i < 3000; i++)
						{
							for (int j = 0; j < 4000 * 3; j = j + 3)
							{
								cameraImage.at<cv::Vec3b>(i, j / 3)[0] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j]);
								cameraImage.at<cv::Vec3b>(i, j / 3)[1] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 1]);
								cameraImage.at<cv::Vec3b>(i, j / 3)[2] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 2]);
							}
						}
						cv::imwrite("cameraImage.bmp", cameraImage);
						num_loop += 1;

						//第一次迭代，获取屏幕图像的定位
						if (num_loop == 1)
						{
							//-------------------------------------处理原理图像------------------------------
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
						}
						else
						{
							////读取相机图像
							//cv::Mat cameraImage;
							//cameraImage = cv::imread("cameraImage.bmp");
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

							//移动标定的中心点
							cv::Point cli_center = cv::Point(2047.22, 1481.77);
							cli_center.x += 2500 - camera_Location.center.x;
							cli_center.y += 2000 - camera_Location.center.y;

							//移动拟合图像
							camera_Location.movedata();

							//获取拟合图像
							camera_Location.getFitdata();

							//相机图像矫正完之后，是相机图像的原始图像了，通过相机图像的拟合图像进行定位，此时拟合坐标已经进行了赋值
							cv::Mat cameraImage_location(4000, 5000, CV_8SC3, cv::Scalar::all(0));

							//cv::Mat cameraImage_location = cv::imread("cameraImage_location.bmp");
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
							if (num_loop == 2)
							{
								camera_Location.generate_P_image();
								right_up = camera_Location.right_up_location_Eye;
								right_down = camera_Location.right_down_location_Eye;
								left_up = camera_Location.left_up_location_Eye;
								left_down = camera_Location.left_down_location_Eye;
							}
							else
							{
								camera_Location.right_up_location_Eye = right_up;
								camera_Location.right_down_location_Eye = right_down;
								camera_Location.left_up_location_Eye = left_up;
								camera_Location.left_down_location_Eye = left_down;
							}



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


							cv::imshow("img", reOriImage);
							cv::waitKey(10);

						}

						cout << "Save a image!!" << endl;
					}

					delete lBufferRGB32;
					lBufferRGB32 = NULL;
				}
			}
			lStream.QueueBuffer(lPtr);
		}

	}

	// Release buffers
	cout << "Releasing buffers" << endl;
	delete[]lBuffers;

	// Now close the stream. Also optionnal but nice to have
	cout << "Closing stream" << endl;
	lStream.Close();

	// Finally disconnect the device. Optional, still nice to have
	cout << "Disconnecting device" << endl;
	lDevice.Disconnect();

	return true;
}

int main()
{
	// PvPipeline used to acquire images from a device
	cout << "1. PvStream sample - image acquisition from a device" << endl << endl;
	// 建立显示窗口
	cv::namedWindow("img", CV_WINDOW_NORMAL);
	cv::setWindowProperty("img", CV_WND_PROP_FULLSCREEN, CV_WINDOW_FULLSCREEN);//设置为全屏

	cv::Mat oriImage, oritest;
	oriImage = cv::imread("ori.bmp");
	oriImage.copyTo(oritest);
	cv::imshow("img", oritest);
	cv::waitKey(10);

	AcquireImages(oriImage);

	cout << endl;
	cout << "<press a key to exit>" << endl;
	PvWaitForKeyPress();

	return 0;
}





