#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);
  
  /**
   * Prediction the next state given that current state is col
   * and time interval is delta_t
   * @param col vector to predict next state for
   * @param delta_t Time between k and k+1 in s
   */
  Eigen::VectorXd PredictionForVector(Eigen::VectorXd & col, double delta_t);
  
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Calculate Lidar measurement for state col
   * @param col state vector calculate measurement from
   */
  Eigen::VectorXd MeasurementLidar(Eigen::VectorXd& col);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Calculate radar measurement for state col
   * @param col state vector calculate measurement from
   */
  Eigen::VectorXd MeasurementRadar(Eigen::VectorXd& col);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  //Lidar measurement matrix
  Eigen::MatrixXd H_;
  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Process noise matrix assembled
  Eigen::MatrixXd P_cov_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Laser measurement covariance matrix assembled
  Eigen::MatrixXd Lidar_M_cov;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Radar measurement covariance matrix assembled
  Eigen::MatrixXd Radar_M_cov;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;
};

#endif  // UKF_H