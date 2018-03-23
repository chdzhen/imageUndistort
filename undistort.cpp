	cv::Mat intrinsicMatrix = cv::Mat::zeros(3, 3, CV_64F);
	cv::Mat distortionCoeff = cv::Mat::ones(1, 4, CV_64F);  
											
	intrinsicMatrix.at<double>(0, 0) = 1598.577215174323,;
	intrinsicMatrix.at<double>(0, 1) = 0;
	intrinsicMatrix.at<double>(0, 2) = 1984.151980970139;
	intrinsicMatrix.at<double>(1, 0) = 0;
	intrinsicMatrix.at<double>(1, 1) = 1600.861242742926;
	intrinsicMatrix.at<double>(1, 2) = 1542.767385830806;
	intrinsicMatrix.at<double>(2, 0) = 0;
	intrinsicMatrix.at<double>(2, 1) = 0;
	intrinsicMatrix.at<double>(2, 2) = 1;

	distortionCoeff.at<double>(0, 0) = -0.00263399;
	distortionCoeff.at<double>(0, 1) = 0.0112513;
	distortionCoeff.at<double>(0, 2) = -0.0163876;
	distortionCoeff.at<double>(0, 3) = 0.00572424;
	


void undistortPoints(cv::InputArray distorted, cv::OutputArray undistorted, 
cv::InputArray K, cv::InputArray D, cv::InputArray R = cv::noArray())
{
	// will support only 2-channel data now for points
	CV_Assert(distorted.type() == CV_32FC2 || distorted.type() == CV_64FC2);
	undistorted.create(distorted.size(), distorted.type());

	//CV_Assert(P.empty() || P.size() == cv::Size(3, 3) || P.size() == cv::Size(4, 3));
	CV_Assert(R.empty() || R.size() == cv::Size(3, 3) || R.total() * R.channels() == 3);
	CV_Assert(D.total() == 4 && K.size() == cv::Size(3, 3) && (K.depth() == CV_32F || K.depth() == CV_64F));

	//��ȡ�ڲ�
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

	//��ȡ�������
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

	//P��ɶ
	//if (!P.empty())
	//{
	//	cv::Matx33d PP;
	//	P.getMat().colRange(0, 3).convertTo(PP, CV_64F);
	//	RR = PP * RR;
	//}

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


void distortPoints(cv::InputArray undistorted, cv::OutputArray distorted, cv::InputArray K, cv::InputArray D, double alpha)
{
	// will support only 2-channel data now for points
	CV_Assert(undistorted.type() == CV_32FC2 || undistorted.type() == CV_64FC2);
	distorted.create(undistorted.size(), undistorted.type());
	size_t n = undistorted.total();

	CV_Assert(K.size() == cv::Size(3, 3) && (K.type() == CV_32F || K.type() == CV_64F) && D.total() == 4);

	//�ڲ�ֵ
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

	//����ֵ
	cv::Vec4d k = D.depth() == CV_32F ? (cv::Vec4d)*D.getMat().ptr<cv::Vec4f>() : *D.getMat().ptr<cv::Vec4d>();

	const cv::Vec2f* Xf = undistorted.getMat().ptr<cv::Vec2f>();
	const cv::Vec2d* Xd = undistorted.getMat().ptr<cv::Vec2d>();
	cv::Vec2f *xpf = distorted.getMat().ptr<cv::Vec2f>();
	cv::Vec2d *xpd = distorted.getMat().ptr<cv::Vec2d>();

	for (size_t i = 0; i < n; ++i)
	{
		cv::Vec2d x = undistorted.depth() == CV_32F ? (cv::Vec2d)Xf[i] : Xd[i];
		//        |x|
		//r2=[x,y]|y| = x2+y2;
		double r2 = x.dot(x);
		double r = std::sqrt(r2);

		// Angle of the incoming ray:
		double theta = atan(r);

		double theta2 = theta*theta, theta3 = theta2*theta, theta4 = theta2*theta2, theta5 = theta4*theta,
			theta6 = theta3*theta3, theta7 = theta6*theta, theta8 = theta4*theta4, theta9 = theta8*theta;

		//����ϵ��
		double theta_d = theta + k[0] * theta3 + k[1] * theta5 + k[2] * theta7 + k[3] * theta9;

		double inv_r = r > 1e-8 ? 1.0 / r : 1;
		double cdist = r > 1e-8 ? theta_d * inv_r : 1;

		cv::Vec2d xd1 = x * cdist;
		cv::Vec2d xd3(xd1[0] + alpha*xd1[1], xd1[1]);
		cv::Vec2d final_point(xd3[0] * f[0] + c[0], xd3[1] * f[1] + c[1]);

		if (undistorted.depth() == CV_32F)
			xpf[i] = final_point;
		else
			xpd[i] = final_point;
	}
}