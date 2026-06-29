# Multi-Sensor Target Tracking via Centralized Kalman Filtering

## 📌 Project Overview
This module implements a **Centralized Multi-Sensor Data Fusion Architecture** utilizing a linear **Kalman Filter (KF)** to track a maneuvering aircraft moving in a two-dimensional plane. By integrating asynchronous or parallel kinematic measurements from multiple independent sensor nodes (e.g., primary and secondary radars), the system minimizes state estimation uncertainty and mitigates the impact of localized measurement noise.

---

<details>
## 📐 Mathematical & System Modeling

The system is framed as a continuous-discrete state-space model tracking the target's 2D position $(x, y)$ and velocities $(\dot{x}, \dot{y})$:

$$\mathbf{x}_k = \begin{bmatrix} x_k & \dot{x}_k & y_k & \dot{y}_k \end{bmatrix}^T$$

### 1. Process Model (Time Update)
Assuming a constant velocity model perturbed by continuous white noise, the state transition is defined by:
$$\mathbf{x}_{k} = \mathbf{A}\mathbf{x}_{k-1} + \mathbf{w}_{k-1}$$

Where $\mathbf{A}$ is the state transition matrix for a sampling time $\Delta t$, and $\mathbf{w}_k \sim \mathcal{N}(0, \mathbf{Q})$ represents the process noise covariance accounting for unmodeled target maneuvers.

### 2. Measurement Model (Measurement Update)
The centralized fusion node receives position coordinates from independent sensors. The measurement vector for sensor $i$ is modeled as:
$$\mathbf{z}_{k}^{i} = \mathbf{H}^{i}\mathbf{x}_k + \mathbf{v}_{k}^{i}$$

Where $\mathbf{H}$ maps the true state space to the measured coordinate space, and $\mathbf{v}_k^{i} \sim \mathcal{N}(0, \mathbf{R}^i)$ represents the unique sensor measurement noise covariance.
</details>
---

## 💻 Code Architecture & Implementation

The module is written purely in **MATLAB** and is split into two primary scripts:

*   **`target_tracking_kf.m`**: The core execution engine. It handles trajectory simulation, sensor noise generation, and executes the recursive Kalman loops (Time Update equations to predict, followed by Centralized Measurement updates to correct).
*   **`test_tracking_pipeline.m`**: A validation script that acts as an automated pipeline test. It executes the tracking engine across varying noise profiles, calculates Root Mean Squared Error (RMSE), and asserts system convergence.

---

## 📊 Key Results & Performance Analysis
*   **Noise Mitigation:** The centralized Kalman filter successfully dampens heavy Gaussian noise injected into raw sensor streams, producing smooth tracking trajectories.
*   **Fusion Gain:** Combining measurements from multiple sensors yields a significantly lower estimation covariance matrix $\mathbf{P}_k$ compared to single-sensor tracking tracking runs, proving the mathematical validity of data fusion.

> 💡 **Detailed Report:** For a comprehensive derivation of the transition matrices, covariance analysis, and full performance plots, refer to the included [tracking_performance.pdf](./tracking_performance.pdf).
