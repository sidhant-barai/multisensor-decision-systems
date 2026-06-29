# Multi-Sensor Target Tracking via Centralized Kalman Filtering

## 📌 Project Overview
This module implements a **Centralized Multi-Sensor Data Fusion Architecture** utilizing a linear **Kalman Filter (KF)** to track a maneuvering aircraft moving in a two-dimensional plane. By integrating tracking data from multiple independent sensors, the system minimizes state estimation uncertainty and mitigates measurement noise.

---

<details>
<summary><b>📐 View Mathematical & System Modeling (Click to expand)</b></summary>

The system tracks the aircraft's 2D position $(x, y)$ and velocities $(\dot{x}, \dot{y})$ using a continuous-discrete state-space model:

```math
\mathbf{x}_k = \begin{bmatrix} x_k & \dot{x}_k & y_k & \dot{y}_k \end{bmatrix}^T
