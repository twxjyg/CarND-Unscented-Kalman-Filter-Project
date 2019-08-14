#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <uWS/uWS.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "json.hpp"
#include "tools.h"
#include "ukf.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;

UKF *ukf_ptr = nullptr;
void HandleCtrl_C(int s);

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // Create a Kalman Filter instance
  UKF ukf;

  ukf_ptr = &ukf;

  // hanle Ctrl-C
  struct sigaction sigIntHandler;

  sigIntHandler.sa_handler = HandleCtrl_C;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;

  sigaction(SIGINT, &sigIntHandler, NULL);

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  h.onMessage([&ukf, &tools, &estimations, &ground_truth](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                                          uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(string(data));

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object
          string sensor_measurement = j[1]["sensor_measurement"];

          MeasurementPackage meas_package;
          std::istringstream iss(sensor_measurement);

          long long timestamp;

          // reads first element from the current line
          string sensor_type;
          iss >> sensor_type;

          if (sensor_type.compare("L") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float px;
            float py;
            iss >> px;
            iss >> py;
            meas_package.raw_measurements_ << px, py;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
          } else if (sensor_type.compare("R") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float ro;
            float theta;
            float ro_dot;
            iss >> ro;
            iss >> theta;
            iss >> ro_dot;
            meas_package.raw_measurements_ << ro, theta, ro_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
          }

          float x_gt;
          float y_gt;
          float vx_gt;
          float vy_gt;
          iss >> x_gt;
          iss >> y_gt;
          iss >> vx_gt;
          iss >> vy_gt;

          VectorXd gt_values(4);
          gt_values(0) = x_gt;
          gt_values(1) = y_gt;
          gt_values(2) = vx_gt;
          gt_values(3) = vy_gt;
          ground_truth.push_back(gt_values);

          // Call ProcessMeasurement(meas_package) for Kalman filter
          try {
            ukf.ProcessMeasurement(meas_package);
          } catch (const std::exception &e) {
            std::cerr << e.what() << std::endl;
            return;
          }

          // Push the current estimated x,y positon from the Kalman filter's
          //   state vector

          VectorXd estimate(4);

          double p_x = ukf.x_(0);
          double p_y = ukf.x_(1);
          double v = ukf.x_(2);
          double yaw = ukf.x_(3);

          double v1 = cos(yaw) * v;
          double v2 = sin(yaw) * v;

          estimate(0) = p_x;
          estimate(1) = p_y;
          estimate(2) = v1;
          estimate(3) = v2;

          estimations.push_back(estimate);

          VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);

          json msgJson;
          msgJson["estimate_x"] = p_x;
          msgJson["estimate_y"] = p_y;
          msgJson["rmse_x"] = RMSE(0);
          msgJson["rmse_y"] = RMSE(1);
          msgJson["rmse_vx"] = RMSE(2);
          msgJson["rmse_vy"] = RMSE(3);
          auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }  // end "telemetry" if

      } else {
        string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket message if
  });  // end h.onMessage

  h.onConnection(
      [&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) { std::cout << "Connected!!!" << std::endl; });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    std::cout << "Disconnected" << std::endl;
    ws.close();
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();

  return 0;
}

void HandleCtrl_C(int s) {
  printf("Caught signal %d\n", s);
  if (ukf_ptr != nullptr) {
    // analyes NIS
    unsigned int target_count = 0;
    std::ofstream nis_file("nis.txt", std::ios::trunc);
    for (auto nis : ukf_ptr->laser_nis_) {
      nis_file << "laser_nis:" << nis << std::endl;
      if (nis > 7.8) {
        target_count++;
      }
    }
    std::cout << "Laser NIS > 7.8 rate:"
              << static_cast<double>(target_count) / static_cast<double>(ukf_ptr->laser_nis_.size()) << std::endl;
    target_count = 0;
    for (auto nis : ukf_ptr->radar_nis_) {
      nis_file << "radar_nis:" << nis << std::endl;
      if (nis > 7.8) {
        target_count++;
      }
    }
    std::cout << "Radar NIS > 7.8 rate:"
              << static_cast<double>(target_count) / static_cast<double>(ukf_ptr->radar_nis_.size()) << std::endl;
    nis_file.close();
  }
  exit(1);
}