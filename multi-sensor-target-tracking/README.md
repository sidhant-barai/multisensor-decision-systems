# 🛰️ Multi-Sensor Flight State Estimation (EKF) Implementation

This workspace contains the implementation scripts, automated testing harness, and sensor bias results for the non-linear 3D flight estimation module.

## 📋 File Inventory

* flight_estimation_ekf.m: The core state estimation engine. It implements an Extended Kalman Filter (EKF) to propagate 3D aircraft dynamics, evaluate non-linear system Jacobians, estimate sensor biases, and isolate data fault anomalies.
* test_tracking_pipeline.m: The automated test and validation harness. It loads the flight data profiles, profiles execution script runtime, and checks for dimension formatting errors or NaN outputs.
* Task4_Bias.png: Graphic output illustrating the tracking performance and stability of the accelerometer and gyroscope bias convergence loops under sensor fault configurations.

## 🚀 How to Run the Implementation

1. Open MATLAB and point your current directory to this specific folder.
2. To verify script dimensions, runtime limits, and basic function execution, run the automated validation test in the command window:

   run('test_tracking_pipeline.m')

3. Ensure that your data file (e.g., dataTask4.mat) is present in your MATLAB path prior to execution.
