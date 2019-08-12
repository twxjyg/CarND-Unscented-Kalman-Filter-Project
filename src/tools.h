#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "domain_types.h"

class Tools {
 public:
  /**
   * A helper method to calculate RMSE.
   */
  static VectorXd CalculateRMSE(const std::vector<VectorXd>& estimations, const std::vector<VectorXd>& ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  static MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
   * A helper method to calcute linear measurement matrix
   */
  static MatrixXd CalculateLinearMeasurementMatrix();

  /**
   * A helper methond to calculate linear State Transform Matrix
   */
  static MatrixXd CalculateLinearStateTransform(const double& det_t);

  /**
   * A helper method to calculate process noise Matrix
   */
  static MatrixXd CalculateProcessNoise(const double& det_t, const double& noise_ax, const double& noise_ay);

  /**
   * A helper method to calculate motion noise Vector
   */
  static VectorXd CalculateMotionNoise(const double& det_t, const double& noise_ax, const double& noise_ay);
  /**
   * A helper method to transform radar measurment to EKF state vector
   */
  static VectorXd TransformRadarMeasurementToState(const VectorXd& radar_measurement);

  /**
   * A helper method to transform laser measurement to EKF state vector
   */
  static VectorXd TransformLaserMeasurementToState(const VectorXd& laser_measurement);

  static MatrixXd GenerateSigmaPoints(const unsigned int& state_dim, const unsigned int& lambda,
                                      const MatrixXd& Psquare_root, const VectorXd& state);

  static MatrixXd PredictSigmaPoints(const unsigned int& state_dim, const double& delta_t, const MatrixXd& Xsig_aug);

  static VectorXd CalculateMeanFromSigmaPoints(const MatrixXd& Xsig_pred, const VectorXd& weights);
  static MatrixXd CalculateCovarianceFromSigmaPoints(const MatrixXd& Xsig_pred, const VectorXd& Xpred,
                                                     const VectorXd& weights);
  static MatrixXd TransformPredictedSigmaPointsToRadarMeasurementSpace(const MatrixXd& Xsig_pred);

  static MatrixXd TransformPredictedSigmaPointsToLaserMeasurementSpace(const MatrixXd& Xsig_pred);
  static MatrixXd CalculateCrossCorrelationMatrix(const unsigned int& state_dim,
                                                  const unsigned int& measurement_dim, const MatrixXd& Zsig,
                                                  const VectorXd& z_pred, const MatrixXd& Xsig_pred, const VectorXd& x,
                                                  const VectorXd& weights);
  /**
   * A helper method to transform state vector to radar measurement space
   */
  static VectorXd TransformStateToRadarMeasurement(const VectorXd& state);

  /**
   * A helper method to transform state vector to laser measurement space
   */
  static VectorXd TransformStateToLaserMeasurement(const VectorXd& state);

  static double NormalizeAngle(const double& phi);
};

#endif  // TOOLS_H_
