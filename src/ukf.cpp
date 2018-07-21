#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);
	x_.fill(0.0);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);


	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1.5;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.5*M_PI;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
	is_initialized_ = false;

	use_laser_ = true;
	use_radar_ = true;

	n_x_ = 5; //state
	n_aug_ = 7; //state of augmentatoion ( add two noise part )
	lambda_ = 3 - n_aug_;

	P_ = MatrixXd::Identity(n_x_, n_x_);


	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1); 
	Xsig_pred_.fill(0.0);
	//after prediction, remove the noise part, so it is 
	//  [n_x,2 * n_aug_ + 1] instead of [n_aug_,2 * n_aug_ + 1]

	weights_ = VectorXd(2*n_aug_+1);//used to calculate mean and covariance
	weights_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	if (!is_initialized_)
	{
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		if ((meas_package.sensor_type_ == MeasurementPackage::LASER) & use_laser_)
		{
			float px, py, v, sai, sai_d;
			px = meas_package.raw_measurements_[0];
			py = meas_package.raw_measurements_[1];
			v = 0;
			sai = 0;
			sai_d = 0;
			x_ << px, py, v, sai, sai_d;
		}
		else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) & use_radar_)
		{
			float ro;
			float theta;
			float ro_dot;

			ro = meas_package.raw_measurements_[0];
			theta = meas_package.raw_measurements_[1];
			ro_dot = meas_package.raw_measurements_[2];

			float px, py,v, sai, sai_d;
			float vx = ro_dot * cos(theta);
			float vy = ro_dot * sin(theta);
			px = ro * cos(theta);
			py = ro * sin(theta);
			v = sqrt(vx*vx + vy*vy);
			sai = 0;
			sai_d = 0;
			x_ << px, py, v, sai, sai_d;
		}

	}
	else
	{
		double delta_t = (meas_package.timestamp_ - time_us_)/ 1000000.0;
		time_us_ = meas_package.timestamp_;
		
		Prediction(delta_t);

		if ((meas_package.sensor_type_ == MeasurementPackage::LASER) & use_laser_)
		{
			UpdateLidar(meas_package);
		}
		else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) & use_radar_)
		{
			UpdateRadar(meas_package);

		}
	}
	
	cout <<endl<< "x_" << endl;
	cout << x_ << endl;
	cout << "============" << endl;
	cout << "P_" << endl;
	cout << P_ << endl;
	cout << "pm finish" << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/
	
	
	/**
	MIAN WORKFLOW:
	part1  generate sigma points
	part2  prediction
	part3  calculate mean and covariance of prediction

	input  : delta_t( used in part2:prediction )
	output : nothing( but modify the x_ and P_ )
	*/


	//=============================== part1.generate sigma points ===========================//
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); //restore sigma points
	Xsig_aug.fill(0.0);

	//mean ( ps: mean of noise == 0 )
	x_aug.head(n_x_) = x_;
	
	//augment of covariance matrix
	MatrixXd Q = MatrixXd(2, 2);
	Q.fill(0.0);
	Q(0, 0) = std_a_*std_a_;
	Q(1, 1) = std_yawdd_*std_yawdd_;

	P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
	P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;
	
	//create sigma points
	Xsig_aug.col(0) = x_aug;
	
	//calculate root of P_aug
	MatrixXd A = MatrixXd(P_.rows(), P_.cols());
	A.fill(0.0);
	A = P_aug.llt().matrixL();
	
	//calculate Xsig_aug :  2 * n_aug + 1 sigma points 
	

	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}
	
	//=============================== part1.generate sigma points ===========================//

	
	//==================================== part2.prediction =================================//
	// according to motion model( in this case CTRV ), predict next time where is the car
	// now Xsig_pred don't include
	float px, py, v, sai, sai_d, mu_a, mu_sai_dd;


	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		px = Xsig_aug(0, i);
		py = Xsig_aug(1, i);
		v = Xsig_aug(2, i);
		sai = Xsig_aug(3, i);
		sai_d = Xsig_aug(4, i);
		mu_a = Xsig_aug(5, i);
		mu_sai_dd = Xsig_aug(6, i);

		if (fabs(sai_d) < 0.001)
		{
			Xsig_pred_(0, i) = px + v*cos(sai)*delta_t + 0.5*delta_t*delta_t*cos(sai)*mu_a;
			Xsig_pred_(1, i) = py + v*sin(sai)*delta_t + 0.5*delta_t*delta_t*sin(sai)*mu_a;
			Xsig_pred_(2, i) = v + 0 + delta_t*mu_a;
			Xsig_pred_(3, i) = sai + 0 + 0.5*delta_t*delta_t*mu_sai_dd;
			Xsig_pred_(4, i) = sai_d + 0 + delta_t*mu_sai_dd;
		}
		else
		{
			Xsig_pred_(0, i) = px + v / sai_d * (sin(sai + sai_d * delta_t) - sin(sai)) + 0.5*delta_t*delta_t*cos(sai)*mu_a;
			Xsig_pred_(1, i) = py + v / sai_d * (-cos(sai + sai_d * delta_t) + cos(sai)) + 0.5*delta_t*delta_t*sin(sai)*mu_a;
			Xsig_pred_(2, i) = v + 0 + delta_t*mu_a;
			Xsig_pred_(3, i) = sai + sai_d*delta_t + 0.5*delta_t*delta_t*mu_sai_dd;
			Xsig_pred_(4, i) = sai_d + 0 + delta_t*mu_sai_dd;
		}

	}

	//==================================== part2.prediction =================================//

	
	//======================= part3.mean and covariance of prediction =======================//
	
	//to calculate mean and covariance, firstly we should know weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);

	for (int i = 1; i < 2*n_aug_+1; i++)
	{
		weights_(i) = 0.5/(lambda_ + n_aug_);
	}
	
	//mean of prediction
	x_.fill(0.0);//every time set zero !!!!
	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	
	//covariance of prediction
	P_.fill(0.0);//every time set zero !!!!
	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		//!!!! careful : make sure angle always in [-Pi , Pi]!!!!
		//reference : soution of "Predicted Mean and Covariance Assignment"
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}

	//======================= part3.mean and covariance of prediction =======================//


	cout << "prediction finish" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
	/**
	MIAN WORKFLOW:
	part1  transform mean and covariance to measurement space ( mean and covariance : from prediction )
	part2  using measurement data to update mean and covariance

	input  : measurement data : radar ( used in part2 )
	output : nothing( but modify the x_ and P_ )
	*/
	cout << "update lidar beginn" << endl;
	int n_z = 2;//lidar date have 2 measurement data
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);
	//======================= part1.transform mean and covariance to measurement space =======================//

	VectorXd z_pred = VectorXd(n_z);//transformed mean
	z_pred.fill(0.0);
	double px, py;
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		px = Xsig_pred_(0, i);
		py = Xsig_pred_(1, i);


		Zsig(0, i) = px;
		Zsig(1, i) = py; 


		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}

	MatrixXd S = MatrixXd(n_z, n_z);//transformed covariance
	S.fill(0.0);
	MatrixXd R = MatrixXd(n_z, n_z);
	R.fill(0.0);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;

		S = S + weights_(i)*(z_diff)*(z_diff).transpose();
	}
	S = S + R;
	//======================= part1.transform mean and covariance to measurement space =======================//


	//====================== part2.using measurement data to update mean and covariance ======================//
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}
	//calculate cross correlation matrix
	MatrixXd K = MatrixXd(n_x_, n_z);
	K.fill(0.0);
	K = Tc*S.inverse();
	//calculate Kalman gain K;

	VectorXd z = VectorXd(n_z);//the measurement data
	z.fill(0.0);
	z = meas_package.raw_measurements_;
	x_ = x_ + K*(z - z_pred);
	


	P_ = P_ - K*S*K.transpose();
	//update state mean and covariance matrix
	//====================== part2.using measurement data to update mean and covariance ======================//

	//====================================== part3.calculate NIS value =======================================//
	double epsilon = (z - z_pred).transpose()*S.inverse()*(z - z_pred);
	cout << "lidar espsilon:"<< epsilon << endl;

	//output to file
	ofstream in;
	in.open("lidar_espsilon.txt", ios::out | ios::app);
	in << epsilon << "\n";
	in.close();
	//====================================== part3.calculate NIS value =======================================//

	cout << "update lidar finish" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/

	/**
	MIAN WORKFLOW:
	part1  transform mean and covariance to measurement space ( mean and covariance : from prediction )
	part2  using measurement data to update mean and covariance

	input  : measurement data : radar ( used in part2 )
	output : nothing( but modify the x_ and P_ )
	*/
	
	cout << "update radar beginn" << endl;
	int n_z = 3;//radar date have 3 measurement data
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);//restore transformed sigam points data
	Zsig.fill(0.0);

	//======================= part1.transform mean and covariance to measurement space =======================//
	VectorXd z_pred = VectorXd(n_z);//transformed mean
	z_pred.fill(0.0);
	double px, py, sai, v;
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		px = Xsig_pred_(0, i);
		py = Xsig_pred_(1, i);
		v = Xsig_pred_(2, i);
		sai = Xsig_pred_(3, i);

		Zsig(0, i) = sqrt(px*px + py*py);
		Zsig(1, i) = atan2(py, px); // must be in [-Pi,Pi]
		Zsig(2, i) = (px*cos(sai)*v + py*sin(sai)*v) / sqrt(px*px + py*py);

		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}
	
	MatrixXd S = MatrixXd(n_z,n_z);//transformed covariance
	S.fill(0.0);
	MatrixXd R = MatrixXd(n_z, n_z);
	R.fill(0.0);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i)*(z_diff)*(z_diff).transpose();
	}
	S = S + R;
	//======================= part1.transform mean and covariance to measurement space =======================//

	
	//====================== part2.using measurement data to update mean and covariance ======================//
	MatrixXd Tc = MatrixXd(n_x_,n_z);
	Tc.fill(0.0);
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
		Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
	}
	//calculate cross correlation matrix
	MatrixXd K = MatrixXd(n_x_, n_z);
	K.fill(0.0);
	K = Tc*S.inverse();
	//calculate Kalman gain K;

	VectorXd z = VectorXd(n_z);//the measurement data
	z.fill(0.0);
	
	z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	x_ = x_ + K*z_diff;

	P_ = P_ - K*S*K.transpose();
	//update state mean and covariance matrix
	//====================== part2.using measurement data to update mean and covariance ======================//
	
	//====================================== part3.calculate NIS value =======================================//
	double epsilon = (z_diff).transpose()*S.inverse()*(z_diff);
	cout << "radar espsilon:" << epsilon << endl;

	//output to file
	ofstream in;
	in.open("radar_espsilon.txt", ios::out | ios::app); 
	in << epsilon << "\n";
	in.close();//¹Ø±ÕÎÄ¼þ

	//====================================== part3.calculate NIS value =======================================//

	
	cout << "update radar finish" << endl;
	
}
