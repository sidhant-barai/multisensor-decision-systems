# 📊 Multi-Objective Optimization and Decision Support Systems

This workspace contains the optimization scripts, evolutionary algorithm parameters, and high-dimensional visualization plots used to design a Pareto-optimal digital control system.

## 📋 File Inventory

* controller_optimization.m: The core decision engine script. It configures and runs an evolutionary algorithm to explore trade-offs between Proportional-Integral (PI) controller gains ($K_p, K_i$).
* parallel_performance_plot.png: High-dimensional visual trade-off representation showing the mapping between optimized controller gains and the system's performance boundaries.
* hypervolume_convergence.png: Metrics convergence graph proving the performance stability and optimization coverage across generations.

## 🚀 How to Run the Implementation

1. Open MATLAB and ensure all evolutionary toolbox dependency scripts are added to your working path directory.
2. Run the main decision framework execution script:

   run('controller_optimization.m')

3. The algorithm will automatically cycle through its generational search phases and render the final system parameter correlations, parallel performance trajectories, and hypervolume stability limits.
