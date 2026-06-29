# 🛰️ Multi-Sensor Target Tracking Implementation

This workspace contains the implementation scripts and analytical reporting for the centralized Kalman filtering target tracking system.

## 📋 File Inventory

* target_tracking_kf.m: The primary execution engine. It generates the simulated 2D aircraft flight path, injects heavy Gaussian sensor noise across multiple independent tracking nodes, and recursively runs the time-prediction and measurement-correction loops.
* test_tracking_pipeline.m: An automated validation and test harness that runs the tracking engine, evaluates filter stability across varying noise thresholds, and logs tracking convergence.
* tracking_performance.pdf: The comprehensive assignment report detailing the kinematic state-space linearization, measurement matrix derivations, innovation sequence tests, and full graphical error analysis plots.

## 🚀 How to Run the Implementation

1. Open MATLAB and set your current working directory to this folder.
2. To run the automated validation suite and verify tracking performance, execute the test harness script in the command window:
   
   run('test_tracking_pipeline.m')

3. To view full coordinate frames, covariance convergence analysis, and localized fault detection systems, open the tracking_performance.pdf file directly.
