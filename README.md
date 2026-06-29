# 📁 Multisensor and Decision Systems Engineering Portfolio

This repository documents advanced implementations of data fusion architectures, kinematic state estimation, and intelligent decision systems. The core focus is on combining information from multiple noisy or asynchronous sensors to achieve robust state tracking, fault detection, and decision-making under uncertainty.

---

## 🛠️ Portfolio Modules

### 1. Multi-Sensor Flight State Estimation via Extended Kalman Filtering

### 📌 Project Overview
This module implements an Advanced Centralized Multi-Sensor Data Fusion Architecture utilizing an Extended Kalman Filter (EKF) to estimate the full 3D kinematic state and sensor biases of an aircraft. By processing parallel measurement streams from high-frequency inertial sensors and GPS data, the system robustly filters measurement noise and isolates sensor fault instances.

---

📐 View Mathematical and System Modeling:

The system tracks the 3D aircraft position, velocity vector, and attitude orientations alongside dynamic gyroscope and accelerometer bias states:

                    x_k = [position_3D; velocity_3D; attitude_3D; biases_6D]

1. Process Model (Time Update)

                               x_k = f(x_k-1, u_k-1) + w_k-1

Where f(x) represents the non-linear rigid-body aircraft flight dynamics, and w_k represents the process noise covariance.

2. Measurement Model (Measurement Update)

                                 z_k = h(x_k) + v_k

Where h(x) maps the state vectors to GPS and air-data measurements (including TAS and Angle of Attack), and v_k represents sensor noise.

---

## 💻 Code Architecture

* flight_estimation_ekf.m: The core data fusion and state estimation engine handling 3D kinematic propagation, non-linear Jacobian evaluations, and recursive EKF corrections.
* test_tracking_pipeline.m: An automated test harness script that loads flight datasets, evaluates dimensions, measures execution runtime, and validates state convergence.

## 📊 Key Results

* Sensor Bias Isolation: Successfully tracked and decoupled structural accelerometer and gyroscope biases from true flight dynamics.
* Robustness to Faults: Maintained flight path state tracking accuracy even during simulated sensor fault anomalies.

💡 Detailed Report: For full performance plots and coordinate transformations, see the included tracking_performance.pdf inside the folder.
