#define _USE_MATH_DEFINES

#include <iostream> 
#include <cmath>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  P_.fill(0.0);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 6;
  
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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  //Lidar measurement matrix
  H_ = MatrixXd(2, 5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;
  P_cov_ = MatrixXd(2, 2);
  P_cov_ << std_a_*std_a_, 0, 
          0, std_yawdd_*std_yawdd_;
  Radar_M_cov = MatrixXd(3, 3);
  Radar_M_cov << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  Lidar_M_cov = MatrixXd(2, 2);
  Lidar_M_cov << std_laspx_*std_laspx_, 0,
                  0, std_laspy_*std_laspy_;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for(int j=1; j < 2*n_aug_+1; j++)
  {
      weights_(j) = 1/(2*(lambda_ + n_aug_));
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_)
  {
    VectorXd meas = meas_package.raw_measurements_;
    if(meas_package.sensor_type_== MeasurementPackage::SensorType::LASER)
    {
      //initiialize 2D position in m
      x_.head(2) = meas;
      //initialize absolute radial velocity in m/s
      x_(3) = 0;
    }
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
      //initiialize 2D position in m
      x_(0) = meas(0)*cos(meas(1));
      x_(1) = meas(0)*sin(meas(1));
      double vx = meas(2)*cos(meas(1));
      double vy = meas(2)*sin(meas(1));
      x_(3) = sqrt(vx*vx + vy*vy);
    }
    //initialize yaw angle
    x_(2) = 0;
    //initialize yaw rate in rad/s
    x_(4) = 0;
    //initialize state covariance matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);
    P_ << 1, 0, 0, 0, 0,
          0, 0.2, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 0.9;                   
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  //calculate time interval in seconds since last measurement update
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  Prediction(delta_t);

  if(meas_package.sensor_type_== MeasurementPackage::SensorType::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  time_us_ = meas_package.timestamp_;
}

VectorXd UKF::PredictionForVector(VectorXd & col, double delta_t)
{
    VectorXd res = VectorXd(5);
    if (col(4) < 0.001 )
    {
       res(0) = col(0) + col(2)*cos(col(3))*delta_t + 0.5*delta_t*delta_t*cos(col(3))*col(5);
        res(1) = col(1) + col(2)*sin(col(3))*delta_t + 0.5*delta_t*delta_t*sin(col(3))*col(5);     
    }
    else
    {
        res(0) = col(0) + col(2)/col(4)*(sin(col(3)+col(4)*delta_t)-sin(col(3))) + 0.5*delta_t*delta_t*cos(col(3))*col(5);
        res(1) = col(1) + col(2)/col(4)*(-cos(col(3)+col(4)*delta_t)+cos(col(3))) + 0.5*delta_t*delta_t*sin(col(3))*col(5);
    }
    res(2) = col(2) + col(5)*delta_t;
    res(3) = col(3) + col(4)*delta_t + 0.5*delta_t*delta_t*col(6);
    res(4) = col(4) + delta_t*col(6);
    return res;  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  //Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  VectorXd y = VectorXd(2);
  y << 0, 0;
  x_aug.tail(2) = y;
  // x_aug(5)=0;
  // x_aug(6)=0;
  //create augmented covariance Matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = P_cov_;
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int j=1; j<=n_aug_; j++)
  {
    Xsig_aug.col(j) = x_aug + sqrt(lambda_+n_aug_)*L.col(j-1);
    Xsig_aug.col(j+n_aug_) = x_aug - sqrt(lambda_+n_aug_)*L.col(j-1);
  }
  //make predictions based on augmented sigma points
  for(int j=0; j<2*n_aug_+1; j++)
  {
    VectorXd col = Xsig_aug.block(0, j, 7, 1);
    Xsig_pred_.col(j) = PredictionForVector(col, delta_t);
  }
  //predict state mean
  x_.fill(0.0);
  for(int j=0; j < 2*n_aug_+1; j++)
  {
      x_ = x_ + weights_(j)*Xsig_pred_.col(j);
  }
  x_(3) = remainder(x_(3), 2*M_PI);
  // predict state covariance matrix
  P_.fill(0.0);
  for(int j=0; j < 2*n_aug_+1; j++)
  {
    VectorXd x_diff = Xsig_pred_.col(j) - x_;
    // angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_ = P_ + weights_(j)*x_diff*x_diff.transpose();  
  }
}

VectorXd UKF::MeasurementRadar(VectorXd& col)
{
    VectorXd res = VectorXd(3);
    res(0) = sqrt(col(0)*col(0) + col(1)*col(1));
    if(col(0) > 0.001) res(1) = atan(col(1)/col(0));
    else{
      if (col(1) > 0) res(1) = M_PI/2;
      else res(1) = -M_PI/2;
    }
    res(2) = (col(0)*cos(col(3))*col(2) + col(1)*sin(col(3))*col(2))/std::max(0.001, res(0));
    return res;
}

VectorXd UKF::MeasurementLidar(VectorXd& col)
{
    VectorXd res = VectorXd(2);
    res(0) = col(0);
    res(1) = col(1);
    return res;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  VectorXd meas = meas_package.raw_measurements_;
  // calculate Kalman gain K;
  MatrixXd Y = (Lidar_M_cov + H_*P_*H_.transpose());
  MatrixXd K = P_*H_.transpose()*Y.inverse();
  // update state mean and covariance matrix
  x_ = x_ + K*(meas-H_*x_);
  //mod out angle by 2*pi to prevent too large values
  x_(3) = remainder(x_(3), 2*M_PI);
  P_ = P_ - K*H_*P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  VectorXd meas = meas_package.raw_measurements_;
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
    // transform sigma points into measurement space
  for(int j=0; j < 2*n_aug_+1; j++)
  {
      VectorXd col = Xsig_pred_.col(j);
      Zsig.col(j) = MeasurementRadar(col);
  }
  // calculate mean predicted measurement
  // std::cout << "z_pred " << z_pred << std::endl;
  for(int j=0; j < 2*n_aug_+1; j++)
  {
      z_pred = z_pred + weights_(j)*Zsig.col(j);
  }
  z_pred(1) = remainder(z_pred(1), 2*M_PI);
  // calculate innovation covariance matrix S
  for(int j=0; j < 2*n_aug_+1; j++)
  {
      S = S + weights_(j)*(Zsig.col(j) - z_pred)*(Zsig.col(j) - z_pred).transpose();
  }
  S = S + Radar_M_cov;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0.0);
  // calculate cross correlation matrix
  for(int j=0; j < 2*n_aug_+1; j++)
  {
      Tc = Tc + weights_(j)*(Xsig_pred_.col(j)-x_)*(Zsig.col(j)-z_pred).transpose();
  }
  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, 3);
  K = Tc*S.inverse();
  // update state mean and covariance matrix
  x_ = x_ + K*(meas-z_pred);
  //mod out angle by 2*pi to prevent too large values
  x_(3) = remainder(x_(3), 2*M_PI);
  P_ = P_ - K*S*K.transpose();
}