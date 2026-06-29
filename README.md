# Multisensor and Decision Systems Engineering Portfolio

This repository documents advanced implementations of data fusion architectures, kinematic state estimation, and intelligent decision systems. The core focus is on combining information from multiple noisy or asynchronous sensors to achieve robust state tracking, fault detection, and decision-making under uncertainty.

---

## Portfolio Modules

### 1. Multi-Sensor Target Tracking via Centralized Kalman Filtering

### Project Overview
This module implements a Centralized Multi-Sensor Data Fusion Architecture utilizing a linear Kalman Filter to track a maneuvering aircraft moving in a two-dimensional plane. By integrating tracking data from multiple independent sensors, the system minimizes state estimation uncertainty and mitigates measurement noise.

---

View Mathematical and System Modeling:

The system tracks the aircraft's 2D position (x, y) and velocities (x_dot, y_dot) using a continuous-discrete state-space model:

x_k = [x_k; x_dot_k; y_k; y_dot_k]

1. Process Model (Time Update)

x_k = A * x_k-1 + w_k-1

Where A is the state transition matrix for a sampling time Delta_t, and w_k represents process noise.

2. Measurement Model (Measurement Update)

z_k^i = H^i * x_k + v_k^i

Where H maps the true state space to the measurements, and v_k^i represents individual sensor noise.

---

## Code Architecture

* target_tracking_kf.m: The core data fusion engine handling trajectory simulation, noise injection, and the recursive Kalman prediction/correction loops.
* test_tracking_pipeline.m: An automated validation script that runs the tracking engine across varying noise profiles to verify convergence and calculate RMSE metrics.

## Key Results

* Fusion Gain: Combining parallel sensor streams significantly lowered the estimation error covariance compared to single-sensor tracking runs.
* Noise Dampening: Effectively filtered high-frequency Gaussian noise to recover smooth, accurate trajectories.

Detailed Report: For full performance plots and coordinate transformations, see the included tracking_performance.pdf inside the folder.
