#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "ukf_exception.h"

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
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  /**
   * accroding to the plotting result under my repo,
   * all three ploting:(NIS_30_30.png, NIS_10_10.png, NIS_9_3.png) are showing that
   * my noise parameter is big enough to describe the measurement noise, but
   * the NIS_9_3 gives me a very good RMSE result, so I think 9 and 3 is enough
   * for now.
   */
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;

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
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  std::cout << meas_package.ToString() << std::endl;
  if (!is_initialized_) {
    x_ = VectorXd::Zero(n_x_);
    P_ = MatrixXd::Identity(n_x_, n_x_) * 1.0;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      x_ = Tools::TransformRadarMeasurementToState(meas_package.raw_measurements_);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ = Tools::TransformLaserMeasurementToState(meas_package.raw_measurements_);
    } else {
      throw UKFException("Unknown measurement type");
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    std::cout << "Initialized X:" << x_.transpose() << std::endl;
    return;
  }

  const double& delta_t = static_cast<double>(meas_package.timestamp_ - time_us_) / 1000000.0;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLaser(meas_package);
    time_us_ = meas_package.timestamp_;
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
    time_us_ = meas_package.timestamp_;
  } else {
    throw UKFException("Unknown measurement type");
  }
  std::cout << "state:" << x_.transpose() << std::endl;
  std::cout << "Cov:" << std::endl << P_ << std::endl;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  // Augmente state vector and state Covariance matrix
  MatrixXd Paug = MatrixXd::Zero(n_aug_, n_aug_);
  Paug.topLeftCorner(n_x_, n_x_) = P_;
  Paug(5, 5) = std::pow(std_a_, 2);
  Paug(6, 6) = std::pow(std_yawdd_, 2);
  lambda_ = 3 - n_aug_;
  MatrixXd Paug_square_toot = Paug.llt().matrixL();
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  // Generate augmented sigma points matrix
  MatrixXd Xsig_aug = Tools::GenerateSigmaPoints(n_aug_, lambda_, Paug_square_toot, x_aug);
  // Predict augmented sigma points
  Xsig_pred_ = Tools::PredictSigmaPoints(n_aug_, delta_t, Xsig_aug);
  // Predict mean
  weights_ = VectorXd::Zero(2 * n_aug_ + 1);
  // 1. set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
  // 2. predict x mean
  x_ = Tools::CalculateMeanFromSigmaPoints(Xsig_pred_, weights_);
  // 3. predict x covariance matrix
  P_ = Tools::CalculateCovarianceFromSigmaPoints(Xsig_pred_, x_, weights_);
  // std::cout << color::red << "pred x:" << x_.transpose() << color::reset << std::endl;
  // std::cout << color::red << "pred P:" << std::endl << P_.transpose() << color::reset << std::endl;
}

void UKF::Update(MeasurementPackage meas_package, const PredictedMeasurementSigmaPointsFunction& GetZsig,
                 const MeasurementNoiseCovarianceFunction& GetR, const NISParamCallback& NISParamCallback) {
  MatrixXd Zsig = GetZsig();
  // std::cout << color::red << "pred Zsig:" << std::endl << Zsig << color::reset << std::endl;
  VectorXd z_pred = Tools::CalculateMeanFromSigmaPoints(Zsig, weights_);
  // std::cout << color::red << "pred z:" << z_pred.transpose() << color::reset << std::endl;
  MatrixXd S = Tools::CalculateCovarianceFromSigmaPoints(Zsig, z_pred, weights_) + GetR();
  if (NISParamCallback != nullptr) {
    NISParamCallback(z_pred, meas_package.raw_measurements_, S);
  }
  // std::cout << color::red << "S:" << std::endl << S << color::reset << std::endl;
  MatrixXd Tc = Tools::CalculateCrossCorrelationMatrix(n_x_, z_pred.rows(), Zsig, z_pred, Xsig_pred_, x_, weights_);
  // std::cout << color::magenta << "Tc:" << std::endl << Tc << color::reset << std::endl;
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
void UKF::UpdateLaser(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  Update(meas_package, [&]() { return Tools::TransformPredictedSigmaPointsToLaserMeasurementSpace(Xsig_pred_); },
         [&]() { return Tools::MakeLaserNoiseMatrix(std_laspx_, std_laspy_); },
         [&](const Eigen::VectorXd& z_pred, const Eigen::VectorXd& z, const Eigen::MatrixXd& S) {
           const double nis = Tools::CalculateNIS(z_pred, z, S);
           std::cout << color::green << "Laser NIS:" << nis << color::reset << std::endl;
           laser_nis_.push_back(nis);
         });
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  Update(meas_package, [&]() { return Tools::TransformPredictedSigmaPointsToRadarMeasurementSpace(Xsig_pred_); },
         [&]() { return Tools::MakeRadarNoiseMatrix(std_radr_, std_radphi_, std_radrd_); },
         [&](const Eigen::VectorXd& z_pred, const Eigen::VectorXd& z, const Eigen::MatrixXd& S) {
           const double nis = Tools::CalculateNIS(z_pred, z, S);
           std::cout << color::green << "Radar NIS:" << nis << color::reset << std::endl;
           radar_nis_.push_back(nis);
         });
}