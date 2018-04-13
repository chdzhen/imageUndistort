// *****************************************************************************
//
//     Copyright (c) 2010, Pleora Technologies Inc., All rights reserved.
//
// *****************************************************************************

//
// To receive images using PvStream directly
//



#include <PvSampleUtils.h>
#include <PvBuffer.h>
#include <PvPixelType.h>
#include <PvBufferWriter.h>
#include <PvString.h>
#include <string>
#include <vector>
#include<io.h>
#include <fstream>

#include"LOCATION1.h"


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



	vector<cv::Point> rowcooder;
	vector<cv::Point> colcooder;
	cv::Point center;

	vector<vector<double>> right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori;
	vector<vector<double>>  right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye;
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
						//相机图像。
						//cameraImage = cv::Mat(3000, 4000, CV_8UC3, cv::Scalar(0, 0, 0));
						for (int i = 0; i < 3000; i++)
						{
							for (int j = 0; j < 4000 * 3; j = j + 3)
							{
								/*cv::Vec3b &cccc = cameraImage.at<cv::Vec3b>(i, j / 3);
								cccc[0] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j]);
								cccc[1] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 1]);
								cccc[2] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 2]);
*/

								cameraImage.at<cv::Vec3b>(i, j/3)[0] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j]);
								cameraImage.at<cv::Vec3b>(i, j / 3)[1] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 1]);
								cameraImage.at<cv::Vec3b>(i, j / 3)[2] = cv::saturate_cast<uchar>(lBuffrtPointr[i * 4000 * 3 + j + 2]);

							}
						}
						cv::imwrite("cameraImage.bmp", cameraImage); 
						cout << "camera have data!!!" << endl;
						/*delete[]lBuffrtPointr;*/
						//lBuffrtPointr = NULL;

						//这里开始获取到cameraimage 
						//locationAndcor(oriImage, cameraImage);
						num_loop += 1;
						if (num_loop == 1)
						{
							vector<cv::Point> othercooder;
							vector<cv::Point> Fit_othercooder;
							vector<vector<double>> right_down_location_fit;
							vector<vector<double>> left_up_location_fit;
							vector<vector<double>> left_down_location_fit;
							vector<vector<double>> right_up_location_fit;
							//原始图像的定位
							ori_fit(oriImage, othercooder, rowcooder, colcooder, center, Fit_othercooder, right_down_location_fit, left_up_location_fit, left_down_location_fit, right_up_location_fit);


							vector<double> RotatecolMargin_JZ;
							vector<double> RotaterowMargin_JZ;
							cv::Point cameracenter_JZ;
							vector<vector<double>> right_down_location, left_up_location, left_down_location, right_up_location;

							//相机图像定位
							camera_location(cameraImage, RotatecolMargin_JZ, RotaterowMargin_JZ, cameracenter_JZ,
								right_down_location, left_up_location, left_down_location, right_up_location);


							//角点对应
							corner_cor(Fit_othercooder, othercooder, RotatecolMargin_JZ, RotaterowMargin_JZ, cameracenter_JZ,
								right_down_location_fit, left_up_location_fit, left_down_location_fit, right_up_location_fit,
								right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori,
								right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye);

							//角点的位置信息
							write_ori_location(right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori);


							//移动点、重新生成oriimage
							reori_show(oriImage, rowcooder, colcooder, center,
								right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori,
								right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye,
								right_down_location, left_up_location, left_down_location, right_up_location);

							cv::destroyAllWindows();

							// 建立显示窗口
							cv::namedWindow("img", CV_WINDOW_NORMAL);
							cv::setWindowProperty("img", CV_WND_PROP_FULLSCREEN, CV_WINDOW_FULLSCREEN);//设置为全屏

							cv::Mat reOriImage;
							reOriImage = cv::imread("reOriImage.jpg");
							cv::imshow("img", reOriImage);
							cv::waitKey(10);
							

						}
						else
						{
							
							vector<double> RotatecolMargin_JZ;
							vector<double> RotaterowMargin_JZ;
							cv::Point cameracenter_JZ;
							vector<vector<double>> right_down_location, left_up_location, left_down_location, right_up_location;

							//相机图像定位
							camera_location(cameraImage, RotatecolMargin_JZ, RotaterowMargin_JZ, cameracenter_JZ,
								right_down_location, left_up_location, left_down_location, right_up_location);

							//角点的位置信息
							write_ori_location(right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori);

							//移动点、重新生成oriimage
							reori_show(oriImage, rowcooder, colcooder, center,
								right_down_location_ori, left_up_location_ori, left_down_location_ori, right_up_location_ori,
								right_down_location_Eye, left_up_location_Eye, left_down_location_Eye, right_up_location_Eye,
								right_down_location, left_up_location, left_down_location, right_up_location);

							cv::Mat reOriImage;
							reOriImage = cv::imread("reOriImage.jpg");
							if (reOriImage.data == NULL)
							{
								cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
							}

							cv::destroyAllWindows();

							// 建立显示窗口
							cv::namedWindow("img", CV_WINDOW_NORMAL);
							cv::setWindowProperty("img", CV_WND_PROP_FULLSCREEN, CV_WINDOW_FULLSCREEN);//设置为全屏

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