#include "tools.h"
#include <iostream>

VectorXd Tools::CalculateRMSE(const std::vector<VectorXd>& estimations, const std::vector<VectorXd>& ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.empty()) {
    return rmse;
  }
  if (estimations.size() != ground_truth.size()) {
    return rmse;
  }

  VectorXd error2_sum(4);
  error2_sum << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd error = estimations[i] - ground_truth[i];
    VectorXd error2 = error.array() * error.array();
    error2_sum += error2;
  }

  error2_sum = error2_sum / estimations.size();

  rmse = error2_sum.array().sqrt();
  return rmse;
}

double Tools::CalculateNIS(const Eigen::VectorXd& z_pred, const Eigen::VectorXd& z, const Eigen::MatrixXd& S) {
  return (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}

VectorXd Tools::TransformRadarMeasurementToState(const VectorXd& radar_measurement) {
  const auto& ro = radar_measurement(0);
  const auto& theta = radar_measurement(1);
  const auto& ro_dot = radar_measurement(2);
  const std::vector<double> data = {ro * std::sin(theta), ro * std::cos(theta), ro_dot, theta, 0.0};
  VectorXd state = VectorXd::Zero(data.size());
  for (unsigned int i = 0; i < data.size(); i++) {
    state(i) = data[i];
  }
  return state;
}

VectorXd Tools::TransformLaserMeasurementToState(const VectorXd& laser_measurement) {
  const std::vector<double>& data = {laser_measurement(0), laser_measurement(1), 0.0, 0.0, 0.0};
  VectorXd state = VectorXd::Zero(data.size());
  for (unsigned int i = 0; i < data.size(); i++) {
    state(i) = data[i];
  }
  return state;
}

MatrixXd Tools::GenerateSigmaPoints(const unsigned int& state_dim, const unsigned int& lambda,
                                    const MatrixXd& Psquare_root, const VectorXd& state) {
  MatrixXd Xsig = MatrixXd::Zero(state_dim, 2 * state_dim + 1);
  // set first column of sigma point matrix
  Xsig.col(0) = state;

  // set remaining sigma points
  for (int i = 0; i < state_dim; ++i) {
    Xsig.col(i + 1) = state + std::sqrt(lambda + state_dim) * Psquare_root.col(i);
    Xsig.col(i + 1 + state_dim) = state - std::sqrt(lambda + state_dim) * Psquare_root.col(i);
  }
  return Xsig;
}

MatrixXd Tools::PredictSigmaPoints(const unsigned int& state_aug_dim, const double& delta_t, const MatrixXd& Xsig_aug) {
  MatrixXd Xsig_pred = MatrixXd(state_aug_dim - 2, 2 * state_aug_dim + 1);
  for (unsigned int i = 0; i < 2 * state_aug_dim + 1; i++) {
    MatrixXd det_x = MatrixXd::Zero(state_aug_dim - 2, 1);
    double v = Xsig_aug.col(i)(2);
    double yaw = Xsig_aug.col(i)(3);
    double yaw_d = Xsig_aug.col(i)(4);
    double Va = Xsig_aug.col(i)(5);
    double Vyaw_dd = Xsig_aug.col(i)(6);
    if (std::abs(yaw_d) < 0.00001) {
      det_x(0) = v * std::cos(yaw) * delta_t;
      det_x(1) = v * std::sin(yaw) * delta_t;
      det_x(2) = 0;
      det_x(3) = yaw_d * delta_t;
      det_x(4) = 0;
    } else {
      det_x(0) = (v / yaw_d) * (std::sin(yaw + yaw_d * delta_t) - std::sin(yaw));
      det_x(1) = (v / yaw_d) * (-std::cos(yaw + yaw_d * delta_t) + std::cos(yaw));
      det_x(2) = 0;
      det_x(3) = yaw_d * delta_t;
      det_x(4) = 0;
    }
    MatrixXd V = MatrixXd::Zero(state_aug_dim - 2, 1);
    V(0) = 0.5 * std::pow(delta_t, 2) * std::cos(yaw) * Va;
    V(1) = 0.5 * std::pow(delta_t, 2) * std::sin(yaw) * Va;
    V(2) = delta_t * Va;
    V(3) = 0.5 * std::pow(delta_t, 2) * Vyaw_dd;
    V(4) = delta_t * Vyaw_dd;
    Xsig_pred.col(i) = Xsig_aug.col(i).head(state_aug_dim - 2) + det_x + V;
  }
  return Xsig_pred;
}

VectorXd Tools::CalculateMeanFromSigmaPoints(const MatrixXd& Xsig_pred, const VectorXd& weights) {
  VectorXd x_pred = VectorXd::Zero(Xsig_pred.rows());
  for (int i = 0; i < Xsig_pred.cols(); ++i) {  // iterate over sigma points
    x_pred = x_pred + weights(i) * Xsig_pred.col(i);
  }
  return x_pred;
}

MatrixXd Tools::CalculateCovarianceFromSigmaPoints(const MatrixXd& Xsig_pred, const VectorXd& Xpred,
                                                   const VectorXd& weights) {
  MatrixXd Ppred = MatrixXd::Zero(Xsig_pred.rows(), Xsig_pred.rows());
  for (int i = 0; i < Xsig_pred.cols(); ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - Xpred;
    // angle normalization
    if (x_diff.rows() == 5) {
      x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    } else if (x_diff.rows() == 3) {
      x_diff(1) = Tools::NormalizeAngle(x_diff(1));
    } else {
      // nothing
    }

    Ppred = Ppred + weights(i) * x_diff * x_diff.transpose();
  }
  return Ppred;
}

MatrixXd Tools::TransformPredictedSigmaPointsToRadarMeasurementSpace(const MatrixXd& Xsig_pred) {
  MatrixXd Zsig = MatrixXd::Zero(3, Xsig_pred.cols());
  for (int i = 0; i < Xsig_pred.cols(); ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double v1 = std::cos(yaw) * v;
    double v2 = std::sin(yaw) * v;

    // measurement model
    Zsig(0, i) = std::sqrt(p_x * p_x + p_y * p_y);                          // r
    Zsig(1, i) = std::atan2(p_y, p_x);                                      // phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / std::sqrt(p_x * p_x + p_y * p_y);  // r_dot
  }
  return Zsig;
}

MatrixXd Tools::TransformPredictedSigmaPointsToLaserMeasurementSpace(const MatrixXd& Xsig_pred) {
  MatrixXd Zsig = MatrixXd::Zero(2, Xsig_pred.cols());
  for (int i = 0; i < Xsig_pred.cols(); ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);

    // measurement model
    Zsig(0, i) = p_x;  // px
    Zsig(1, i) = p_y;  // py
  }
  return Zsig;
}

MatrixXd Tools::CalculateCrossCorrelationMatrix(const unsigned int& state_dim, const unsigned int& measurement_dim,
                                                const MatrixXd& Zsig, const VectorXd& z_pred, const MatrixXd& Xsig_pred,
                                                const VectorXd& x, const VectorXd& weights) {
  assert(Xsig_pred.cols() == Zsig.cols());
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(state_dim, measurement_dim);

  // calculate cross correlation matrix
  for (int i = 0; i < Xsig_pred.cols(); ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    if (z_diff.rows() == 3) {
      z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  return Tc;
}

MatrixXd Tools::MakeRadarNoiseMatrix(const double& std_radius, const double& std_angle, const double& std_radius_d) {
  MatrixXd R = MatrixXd::Zero(3, 3);
  R(0, 0) = std::pow(std_radius, 2);
  R(1, 1) = std::pow(std_angle, 2);
  R(2, 2) = std::pow(std_radius_d, 2);
  return R;
}

MatrixXd Tools::MakeLaserNoiseMatrix(const double& std_px, const double& std_py) {
  MatrixXd R = MatrixXd::Zero(2, 2);
  R(0, 0) = std::pow(std_px, 2);
  R(1, 1) = std::pow(std_py, 2);
  return R;
}

double Tools::NormalizeAngle(const double& phi) { return std::atan2(std::sin(phi), std::cos(phi)); }