#include "ukf.h"
#include "Eigen/Dense"
#include "math.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  Xsig_pred_ = MatrixXd(5,15);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */
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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  time_us_ = 0;
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

	x_ << 0,0,0,0,0;
	P_ = MatrixXd::Identity(5,5);
	P_(3,3) = 0.0225;
	P_(4,4) = 0.0225;

	n_x = 5;
	n_z = 3; // radar measurement dimension
	n_aug = 7;
	lambda = 3 - n_aug;

	weights = VectorXd(n_aug*2 + 1);

	weights(0) = lambda / (lambda + n_aug);
	for (int i=1; i<2*n_aug+1; i++)
	{
	    weights(i) = 1/(2*(lambda+n_aug));
	}

is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
	double delta_t;
	delta_t = static_cast<double>(meas_package.timestamp_ - time_us_)/1000000;
	// cout << "Prediction done "<<endl;
	// if (is_initialized_==false)
	// {
	// 	if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	// 	{
	// 		x_(0) = meas_package.raw_measurements_(0);
	// 		x_(1) = meas_package.raw_measurements_(1);
	// 	}
	// 	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	// 	{
	// 		x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
	// 		x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
	// 		x_(2) = meas_package.raw_measurements_(2);
	// 	}
	// 	time_us_ = meas_package.timestamp_;
	// 	is_initialized_ = true;
	// 	return;
	// }

	if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double rho_d = meas_package.raw_measurements_(2);
            double p_x = rho * cos(phi);
            double p_y = rho * sin(phi);
            // To note: this is inaccurate value, aim to initialize velocity which's magnitude/order is close to real value
            double vx = rho_d * cos(phi);
            double vy = rho_d * sin(phi);
            double v = sqrt(vx * vx + vy * vy);
            x_ << p_x, p_y, v,0, 0;
           
        } else {
            x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0., 0, 0;
            
        }
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
  	}	

	Prediction(delta_t);
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
	{
		UpdateLidar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
	{
		UpdateRadar(meas_package);		
	}
	time_us_ = meas_package.timestamp_;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
	// int n_x = 5, n_aug = 7;
	// int lambda = 3 - n_aug;

	MatrixXd Xsig_aug(n_aug,2*n_aug+1);

	VectorXd x_aug(n_aug);
	x_aug.fill(0.0);
	x_aug.head(n_x) = x_;

	MatrixXd P_aug(n_aug,n_aug);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x,n_x) = P_;
	P_aug(n_x,n_x) = std_a_ * std_a_;
	P_aug(n_x+1,n_x+1) = std_yawdd_ * std_yawdd_;
	// MatrixXd Q(2,2);
	// Q << std_a_*std_a_, 0,
	// 	0 , std_yawdd_*std_yawdd_;
	// P_aug.bottomRightCorner(2,2) = Q;

	// cout << P_aug << endl;

	MatrixXd L = P_aug.llt().matrixL();

	Xsig_aug.col(0) = x_aug;
	for (int i=0; i<n_aug; i++)
	{
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda + n_aug) * L.col(i);
		Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i);
	}

	double px, py, v, yaw, yaw_rate, nu_a, nu_yaw;
  
	for (int i=0; i<2*n_aug+1; i++)
	{
		px = Xsig_aug(0,i); py = Xsig_aug(1,i); v = Xsig_aug(2,i);
		yaw = Xsig_aug(3,i); yaw_rate = Xsig_aug(4,i);
		nu_a = Xsig_aug(5,i); nu_yaw = Xsig_aug(6,i);
	  
		if (yaw_rate == 0)
		{
			Xsig_pred_(0,i) = Xsig_aug(0,i) + v*cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
			Xsig_pred_(1,i) = Xsig_aug(1,i) + v*sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
			Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t * nu_a;
			Xsig_pred_(3,i) = Xsig_aug(3,i) + 0.5*delta_t*delta_t*nu_yaw;
			Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t*nu_yaw;
		}

		else
		{
			Xsig_pred_(0,i) = Xsig_aug(0,i) + v*(sin(yaw + yaw_rate*delta_t) - sin(yaw))/yaw_rate + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
			Xsig_pred_(1,i) = Xsig_aug(1,i) + v*(-cos(yaw+yaw_rate*delta_t) + cos(yaw))/yaw_rate + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
			Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t * nu_a;
			Xsig_pred_(3,i) = Xsig_aug(3,i) + yaw_rate*delta_t + 0.5*delta_t*delta_t*nu_yaw;
			Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t*nu_yaw;
		}
	}

    VectorXd temp(n_x);
    temp.fill(0.0);
    for (int i=0; i<2*n_aug+1; i++)
    {
        temp = temp + weights(i) * Xsig_pred_.col(i);
    }
  // predict state covariance matrix
    MatrixXd Ptemp(n_x,n_x);
    Ptemp.fill(0.0);
    
    for (int i=0; i<2*n_aug+1; i++)
    {
        // Ptemp = (Xsig_pred_.col(i)-x);
        // Ptemp = weights(i) * Ptemp * Ptemp.transpose(); 
        // P = P + Ptemp;
        Ptemp = Ptemp + weights(i) * (Xsig_pred_.col(i)-temp) * (Xsig_pred_.col(i)-temp).transpose();
    }

    x_ = temp;
    P_ = Ptemp;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
	
	MatrixXd R(2,2);

	VectorXd y(2), z(2);
	MatrixXd H(2,n_x), S(2,2), K(n_x,2), I(n_x,n_x);

	I = MatrixXd::Identity(n_x,n_x);
	R << std_laspx_ * std_laspx_, 0, 
		0 , std_laspy_ * std_laspy_;
	H << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;

	z = meas_package.raw_measurements_;
	y = z - H * x_;
	S = H * P_ * H.transpose() + R;
	K = P_ * H.transpose() * S.inverse();

	x_ = x_ + K * y;
	P_ = (I - K * H) * P_;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
	// int n_x=5, n_z=3, lambda=3-n_x;

	VectorXd z(n_z);
	z = meas_package.raw_measurements_;

  	MatrixXd Zsig(n_z,2*n_aug+1);
	double px, py, v, yaw, yaw_rate;
    for (int i=0; i<2*n_aug+1; i++)
    {
        px = Xsig_pred_(0,i); py = Xsig_pred_(1,i); v = Xsig_pred_(2,i);
        yaw = Xsig_pred_(3,i); yaw_rate = Xsig_pred_(4,i);
        
        Zsig(0,i) = sqrt(px*px + py*py);
        Zsig(1,i) = atan2(py,px);
        Zsig(2,i) = (px*v*cos(yaw) + py*v*sin(yaw))/Zsig(0,i);
    }
  
    // VectorXd weights(2*n_x+1);
    // weights(0) = lambda/(lambda+n_x);
    // for (int i=1; i<2*n_x+1; i++)
    // {
    // 	weights(i) = 1/(2*(lambda+n_x));
    // }

    VectorXd z_pred(n_z);
	z_pred.fill(0.0);
    for (int i=0;i<2*n_aug+1; i++)
    {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    // if (z_pred(1)>M_PI)  z_pred(1) -= 2*M_PI;
    // else if (z_pred(1)<-M_PI)  z_pred(1) += 2*M_PI;

    VectorXd zdiff(n_z);
	MatrixXd S(n_z,n_z);
    S.fill(0.0);
    for (int i=0; i<2*n_aug+1; i++)
    {
        zdiff = Zsig.col(i) - z_pred;
	    if (zdiff(1)>M_PI) zdiff(1) -= 2*M_PI;
	    if (zdiff(1)<-M_PI) zdiff(1) += 2*M_PI;
        S = S + weights(i) * zdiff * zdiff.transpose() ;
    }
    
    MatrixXd R(n_z,n_z);
    R << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  
    S = S + R;
  
    VectorXd xdiff(n_x);
    MatrixXd Tc(n_x,n_z);
    Tc.fill(0.0);
    for (int i=0; i<2*n_aug+1;i++)
    {
	    xdiff = Xsig_pred_.col(i) - x_;
	    zdiff = Zsig.col(i) - z_pred;
	    if (xdiff(3)>M_PI) xdiff(3) -= 2*M_PI;
	    if (xdiff(3)<-M_PI) xdiff(3) += 2*M_PI;
	    if (zdiff(1)>M_PI) zdiff(1) -= 2*M_PI;
	    if (zdiff(1)<-M_PI) zdiff(1) += 2*M_PI;
    	Tc = Tc + weights(i) * (xdiff) * (zdiff).transpose();
    }

    MatrixXd K(n_x,n_z);
    K = Tc * S.inverse();

    x_ = x_ + K * (z - z_pred);
    P_ = P_ - K * S * K.transpose();
}