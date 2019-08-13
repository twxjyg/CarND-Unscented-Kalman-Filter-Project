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
   * A helper method to calculate NIS(Normalized Innovation Squared)
   */
  static double CalculateNIS(const Eigen::VectorXd& z_pred, const Eigen::VectorXd& z, const Eigen::MatrixXd& S);

  /**
   * A helper method to transform radar measurment to KF state vector
   */
  static VectorXd TransformRadarMeasurementToState(const VectorXd& radar_measurement);

  /**
   * A helper method to transform laser measurement to KF state vector
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

  static MatrixXd CalculateCrossCorrelationMatrix(const unsigned int& state_dim, const unsigned int& measurement_dim,
                                                  const MatrixXd& Zsig, const VectorXd& z_pred,
                                                  const MatrixXd& Xsig_pred, const VectorXd& x,
                                                  const VectorXd& weights);

  static MatrixXd MakeRadarNoiseMatrix(const double& std_radius, const double& std_angle, const double& std_radius_d);

  static MatrixXd MakeLaserNoiseMatrix(const double& std_px, const double& std_py);

  static double NormalizeAngle(const double& phi);
};

#endif  // TOOLS_H_
