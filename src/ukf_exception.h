#include <stdexcept>
#ifndef UKF_EXCEPTION_UKF_EXCEPTION_H_
#define UKF_EXCEPTION_UKF_EXCEPTION_H_
class UKFException : public std::runtime_error {
 public:
  explicit UKFException(const std::string& what) : std::runtime_error(what) {}
};

#endif  // UKF_EXCEPTION_UKF_EXCEPTION_H_