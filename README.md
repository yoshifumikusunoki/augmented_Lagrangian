# Introduction
MATLAB sources for the augmented Lagrangian method

# Execution

Run the following.
```
exp_AL_quad
```
It solve a convex quadratic problem whose parameters are read from `prob_quad_var2_ineq4.mat`.
The result is saved in `result_AL.mat`.
The source of the augmented Lagrangian is `augmented_Lagrangian.m`.

Run the following.
```
plot_prob_quad
```
Then, you see a plot of the quadratic problem and the trajectory obtained by the augmented Lagrangian.
