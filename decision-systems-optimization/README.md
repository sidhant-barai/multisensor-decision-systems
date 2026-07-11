# 📊 Multi-Objective Optimization and Decision Support Systems

This workspace documents an evolutionary computational design engine focused on high-dimensional search optimization and multi-criterion decision analysis (MCDA). The framework systematically resolves engineering trade-offs between transient performance criteria and control energy bounds based on explicit stakeholder goals.

## 📋 File Inventory

* controller_optimization.m: The core decision framework execution script. It initializes search spaces via multi-strategy plans (Full Factorial, Latin Hypercube, and Sobol sequences) and runs a multi-objective evolutionary loop to evaluate digital PI controller gains ($K_p, K_i$).
* parallel_tradeoffs.png: High-dimensional parallel coordinate trajectory plot showing the mapping between optimized controller gains and the system's performance boundaries.
* hypervolume_convergence.png: Metrics convergence graph tracking optimization coverage and front stability over 50 generations.
* criterion_correlations.png: Matrix scatter plot showing performance criteria trade-offs and structural parameter correlations post-optimization.
