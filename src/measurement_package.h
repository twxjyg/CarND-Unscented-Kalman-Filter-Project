#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "domain_types.h"
#include "color.h"

class MeasurementPackage {
 public:
  enum SensorType { LASER, RADAR } sensor_type_;

  long long timestamp_;

  VectorXd raw_measurements_;
  std::string ToString() const {
    std::stringstream ss;
    if (sensor_type_ == LASER) {
      ss << color::red << "laser:(px=" << raw_measurements_(0) << ","
         << "py=" << raw_measurements_(1) << ")" << color::reset;
    } else {
      ss << color::green << "radar:(ro=" << raw_measurements_(0) << ","
         << "theta=" << raw_measurements_(1) << ","
         << "ro_d=" << raw_measurements_(2) << ")" << color::reset;
    }
    return ss.str();
  }
};

#endif  // MEASUREMENT_PACKAGE_H_