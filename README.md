# Numerical Simulation of Heat Conduction using Crank-Nicolson Method
A numerical simulation of heat conduction using the Crank-Nicolson implicit method. The heat equation is discretized in both time and space, and the resulting tridiagonal system is solved using the Thomas algorithm. The code solves the heat equation numerically and compares the results with the analytical solution. The simulation is conducted over a specified number of time steps and spatial steps.

## Requirements

```
gnuplot (plot)
gcc (compile and run)
```

# Description

The heat equation is given by $$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$ where *T* is the temperature, *t* is time, *x* is position, and *α* is the thermal diffusivity. The one-dimensional equation is considered for this simulation.

## Crank-Nicolson Implicit Method

The Crank-Nicolson method is used to discretize the heat equation in
time and space. The discretized equation is given by
$$\frac{T\_{i}^{j+1} - T\_{i}^{j}}{\Delta t} = \frac{\alpha}{2} \left(\frac{T\_{i+1}^{j+1} - 2T\_{i}^{j+1} + T\_{i-1}^{j+1}}{\Delta x^2} + \frac{T\_{i+1}^{j} - 2T\_{i}^{j} + T\_{i-1}^{j}}{\Delta x^2}\right)$$
where: $T\_{i}^{j} \text{ is the temperature at spatial position}$ $x_i \text{ and time step } t_j$, $\Delta t \text{ is the time step size,}$ $\Delta x \text{ is the spatial step size, and}$ $\alpha \text{ is the thermal diffusivity.}$
This equation represents the change in temperature over time due to the diffusion of heat.

## Parameter *r*

In the Crank-Nicolson method, the parameter *r* is defined as:
$$r = \frac{\alpha \Delta t}{\Delta x^2}$$
It represents the ratio of the heat diffusion rate to the spatial step
size squared. The value of *r* is crucial in determining the accuracy of
the numerical solution. Unlike explicit methods, implicit methods like
Crank-Nicolson are unconditionally stable, meaning they can handle
larger time step sizes without leading to instability.

## Tridiagonal System and Thomas Algorithm

The discretized equation results in a tridiagonal system of equations.
This system can be represented as follows:
```math
$$\begin{bmatrix}
b_1 & c_1 & 0 & 0 & \cdots & 0 \\
a_2 & b_2 & c_2 & 0 & \cdots & 0 \\
0 & a_3 & b_3 & c_3 & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & 0 & \cdots & b_n
\end{bmatrix}
\begin{bmatrix}
T_1^{j+1} \\
T_2^{j+1} \\
T_3^{j+1} \\
\vdots \\
T_n^{j+1}
\end{bmatrix}
=
\begin{bmatrix}
d_1 \\
d_2 \\
d_3 \\
\vdots \\
d_n
\end{bmatrix}$$
```
where $a_i = -1$, $b_i = 2 \left(\frac{1}{r} + 1\right)$, $c_i = -1$,
and
$d_i = T_{i-1}^{j} + 2 \left(\frac{1}{r} - 1\right) T_{i}^{j} + T_{i+1}^{j}$.
The values of $T$ are solved using the Thomas algorithm involving
forward elimination and back substitution.

## Coefficients

The coefficients for the tridiagonal system are given by:
```math
$$\begin{aligned}
a_i &= -1 \\
b_i &= 2 \left(\frac{1}{r} + 1\right) \\
c_i &= -1 \\
d_i &= T_{i-1}^{j} + 2 \left(\frac{1}{r} - 1\right) T_{i}^{j} + T_{i+1}^{j}
\end{aligned}$$
```

## Boundary Conditions

With the following boundary conditions: 
```math
$$\begin{aligned}
d_1 &= T_0 & d_n &= T_L \\
b_1 &= 1 & a_n &= 0 \\
c_1 &= 0 & b_n &= 1
\end{aligned}$$
```

Taking the values, we get: 
```math
$$\begin{aligned}
a &= \begin{bmatrix} 0 \\ -1 \\ -1 \\ -1 \end{bmatrix}, &
b &= \begin{bmatrix} 1 \\ 2 \left(\frac{1}{r} + 1\right) \\
2 \left(\frac{1}{r} + 1\right) \\ 2 \left(\frac{1}{r} + 1\right) \\ 1 \end{bmatrix}, &
c &= \begin{bmatrix} -1 \\ -1 \\ -1 \\ 0 \end{bmatrix}, &
d &= \begin{bmatrix} T_{0} \\ T_{1} + 2 \left(\frac{1}{r} - 1\right) T_{2} + T_{3} \\
T_{2} + 2 \left(\frac{1}{r} - 1\right) T_{3} + T_{4} \\ T_{3} + 2 \left(\frac{1}{r} - 1\right) T_{4} + T_{5} \\
T_{L} \end{bmatrix}
\end{aligned}$$
```

The final tridiagonal system matrix multiplication can be represented
as:
```math
$$\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\
-1 & 2 \left(\frac{1}{r} + 1\right) & -1 & 0 & 0 \\
0 & -1 & 2 \left(\frac{1}{r} + 1\right) & -1 & 0 \\
0 & 0 & -1 & 2 \left(\frac{1}{r} + 1\right) & -1 \\
0 & 0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
T_1^{j+1} \\
T_2^{j+1} \\
T_3^{j+1} \\
T_4^{j+1} \\
T_5^{j+1}
\end{bmatrix}
=
\begin{bmatrix}
T_{0} \\
T_{1}^{j} + 2 \left(\frac{1}{r} - 1\right) T_{2}^{j} + T_{3}^{j} \\
T_{2}^{j} + 2 \left(\frac{1}{r} - 1\right) T_{3}^{j} + T_{4}^{j} \\
T_{3}^{j} + 2 \left(\frac{1}{r} - 1\right) T_{4}^{j} + T_{L}^{j} \\
T_{L}
\end{bmatrix}$$
```

Solving this tridiagonal system using the Thomas algorithm gives us the
values of $T_i^{j+1}$, which represent the temperature distribution at
the next time step.
