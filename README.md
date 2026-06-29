# Multi-Sensor Target Tracking via Centralized Kalman Filtering

## 📌 Project Overview
This module implements a **Centralized Multi-Sensor Data Fusion Architecture** utilizing a linear **Kalman Filter (KF)** to track a maneuvering aircraft moving in a two-dimensional plane. By integrating tracking data from multiple independent sensors, the system minimizes state estimation uncertainty and mitigates measurement noise.

---

<details>
<summary><b>📐 View Mathematical & System Modeling (Click to expand)</b></summary>

The system tracks the aircraft's 2D position $(x, y)$ and velocities $(\dot{x}, \dot{y})$ using a continuous-discrete state-space model:

$$\mathbf{x}_k = \begin{bmatrix} x_k & \dot{x}_k & y_k & \dot{y}_k \end{bmatrix}^T$$

### 1. Process Model (Time Update)
$$\mathbf{x}_{k} = \mathbf{A}\mathbf{x}_{k-1} + \mathbf{w}_{k-1}$$
Where $\mathbf{A}$ is the state transition matrix for a sampling time $\Delta t$, and $\mathbf{w}_k \sim \mathcal{N}(0, \mathbf{Q})$ represents process noise.

### 2. Measurement Model (Measurement Update)
$$\mathbf{z}_{k}^{i} = \mathbf{H}^{i}\mathbf{x}_k + \mathbf{v}_{k}^{i}$$
Where $\mathbf{H}$ maps the true state space to the measurements, and $\mathbf{v}_k^{i} \sim \mathcal{N}(0, \mathbf{R}^i)$ represents individual sensor noise.
</details>

---

## 💻 Code Architecture
* **`target_tracking_kf.m`**: The core data fusion engine handling trajectory simulation, noise injection, and the recursive Kalman prediction/correction loops.
* **`test_tracking_pipeline.m`**: An automated validation script that runs the tracking engine across varying noise profiles to verify convergence and calculate RMSE metrics.

## 📊 Key Results
* **Fusion Gain:** Combining parallel sensor streams significantly lowered the estimation error covariance ($\mathbf{P}_k$) compared to single-sensor tracking runs.
* **Noise Dampening:** Effectively filtered high-frequency Gaussian noise to recover smooth, accurate trajectories.

> 💡 **Detailed Report:** For full performance plots and coordinate transformations, see the included [tracking_performance.pdf](./tracking_performance.pdf).
