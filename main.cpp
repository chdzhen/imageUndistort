
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

	//�洢����ͼ������꼰λ����Ϣ
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

						//��һ�ε�������ȡ��Ļͼ��Ķ�λ
						if (num_loop == 1)
						{
							//-------------------------------------����ԭ��ͼ��------------------------------
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
						}
						else
						{
							////��ȡ���ͼ��
							//cv::Mat cameraImage;
							//cameraImage = cv::imread("cameraImage.bmp");
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

							//�ƶ��궨�����ĵ�
							cv::Point cli_center = cv::Point(2047.22, 1481.77);
							cli_center.x += 2500 - camera_Location.center.x;
							cli_center.y += 2000 - camera_Location.center.y;

							//�ƶ����ͼ��
							camera_Location.movedata();

							//��ȡ���ͼ��
							camera_Location.getFitdata();

							//���ͼ�������֮�������ͼ���ԭʼͼ���ˣ�ͨ�����ͼ������ͼ����ж�λ����ʱ��������Ѿ������˸�ֵ
							cv::Mat cameraImage_location(4000, 5000, CV_8SC3, cv::Scalar::all(0));

							//cv::Mat cameraImage_location = cv::imread("cameraImage_location.bmp");
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
	// ������ʾ����
	cv::namedWindow("img", CV_WINDOW_NORMAL);
	cv::setWindowProperty("img", CV_WND_PROP_FULLSCREEN, CV_WINDOW_FULLSCREEN);//����Ϊȫ��

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





