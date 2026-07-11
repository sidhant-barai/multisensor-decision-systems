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

### 2. Multi-Objective Optimization and Decision Support Systems

### 📌 Project Overview
This module implements an evolutionary computational design engine focused on search optimization and multi-criterion decision analysis (MCDA). The system searches a high-dimensional parameter space to identify Pareto-optimal digital controller gains, systematically resolving engineering trade-offs between system performance criteria and control energy bounds based on stakeholder preferences.

---

📐 View Mathematical and Optimization Modeling:

The system maps the vector of decision variables (controller gains) to a multi-dimensional objective function space:

                          f(gains) = [J_1, J_2, J_3, J_4]^T

Where the algorithm minimizes conflicting transient performance indexes and actuation costs simultaneously:

             Minimize:  J_1 (Rise Time),  J_2 (Overshoot),  J_3 (Error),  J_4 (Control Effort)

Subject to preferential boundaries defined by explicit performance goals:

                                   J_i <= Goal_i

---

## 💻 Code Architecture

* 📁 01-multi-sensor-flight-estimation/: Core estimation engine handling 3D kinematic propagation, Jacobian evaluations, and recursive EKF noise filtering.
* 📁 02-decision-systems-optimization/: Evolutionary multi-objective search scripts, trade-off analysis code, and hypervolume metric convergence tracking.
