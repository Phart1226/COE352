# COE 352 Project 2 - Finite Element Solver
### _Preston Hart_

Link to Github repo: [Project 2](https://github.com/Phart1226/COE352/tree/master/Project2)

### Prereques For Running the Code
- pip install numpy, sympy, and matplotlib
- sympy was used for debugging the code to look at symbolic representations of the integrations taking place

## Project Overview
For this Finite Element solver, the heat equation (ut - uxx = f(x,t)) was the focus. In general, the finite element process constructs a mass and stiffness matrix from Lagragian basis functions that represent the connection of nodes and elements in a global space.
This specific project was built using a 1D Galerkin structure, meaning 1D Lagrange basis functions were used in integration to create the mass and stiffness matrices. The specific function f(x,t) = (pi^2 -1)e^(-t)sin(pi*x) and u(x,0) = sin(pi*x)

### Running the Code
After intial startup the user should see a printout like:

```
Welcome to the Finite Element Simulator
Please fill in the following user input:


Number of nodes:
```


For this example and project, 11 nodes, and therefore 10 elements, will be used.

The following is an example of a correctly formatted user input:
```
Welcome to the Finite Element Simulator
Please fill in the following user input:


Number of nodes: 11
Location of last node: 1
Boundary condition 1: 0
Boundary condition 2: 0
Boundary condition 3: 1
Boundary condition 4: 1
Boundary condition 5: 0
Boundary condition 6: 0
Number of timesteps: 600
Symbolic representation for debugging?(True or False) False
FE or BE? FE
```
- Location of last node: This code assumes that the first node starts at location 0 and the integer value entered here will determine the width of each element
- Boundary condition 1 and 2: u(0,t) and u(1,t)
- Boundary condition 3 and 4: Derichlet Boundary conditions applied
- Boundary condition 5 and 6: u(0,0) and u(1,1)
- Number of timesteps: 1/dt, this project is assuming that the final time is 1
- Symbolic representation: This feature turns sympy representations of the matrices that are created on or off. For numercial results to be found, False should always be entered here
- FE or BE? Forward Euler or Backward Euler numerical discretization method

### Creating the Mass and Stiffness Matrices and Force Vector

#### Elemental Matrices/Vector
Mass:
```
[[h/3, h/6], 
 [h/6, h/3]]
```
Stiffness
```
[[1/h, -1/h], 
 [-1/h, 1/h]]
 ```

 Force Vector
 ```
 [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X)], 
 [h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)]
 ```



The mass and stiffness matrices that are created from the above input using a symbolic representation: 

```
[[h/3, h/6, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
 [h/6, 2*h/3, h/6, 0, 0, 0, 0, 0, 0, 0, 0], 
 [0, h/6, 2*h/3, h/6, 0, 0, 0, 0, 0, 0, 0], 
 [0, 0, h/6, 2*h/3, h/6, 0, 0, 0, 0, 0, 0], 
 [0, 0, 0, h/6, 2*h/3, h/6, 0, 0, 0, 0, 0], 
 [0, 0, 0, 0, h/6, 2*h/3, h/6, 0, 0, 0, 0], 
 [0, 0, 0, 0, 0, h/6, 2*h/3, h/6, 0, 0, 0], 
 [0, 0, 0, 0, 0, 0, h/6, 2*h/3, h/6, 0, 0], 
 [0, 0, 0, 0, 0, 0, 0, h/6, 2*h/3, h/6, 0], 
 [0, 0, 0, 0, 0, 0, 0, 0, h/6, 2*h/3, h/6], 
 [0, 0, 0, 0, 0, 0, 0, 0, 0, h/6, h/3]]


 [[1/h, -1/h, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
  [-1/h, 2/h, -1/h, 0, 0, 0, 0, 0, 0, 0, 0], 
  [0, -1/h, 2/h, -1/h, 0, 0, 0, 0, 0, 0, 0], 
  [0, 0, -1/h, 2/h, -1/h, 0, 0, 0, 0, 0, 0], 
  [0, 0, 0, -1/h, 2/h, -1/h, 0, 0, 0, 0, 0], 
  [0, 0, 0, 0, -1/h, 2/h, -1/h, 0, 0, 0, 0], 
  [0, 0, 0, 0, 0, -1/h, 2/h, -1/h, 0, 0, 0], 
  [0, 0, 0, 0, 0, 0, -1/h, 2/h, -1/h, 0, 0], 
  [0, 0, 0, 0, 0, 0, 0, -1/h, 2/h, -1/h, 0], 
  [0, 0, 0, 0, 0, 0, 0, 0, -1/h, 2/h, -1/h], 
  [0, 0, 0, 0, 0, 0, 0, 0, 0, -1/h, 1/h]]
 ```

 Force Vector:
 ```
 [[h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(1/2 - X/2)*exp(-t)*sin(pi*X) + h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)], 
  [h*(-1 + pi**2)*(X/2 + 1/2)*exp(-t)*sin(pi*X)]]
```
However, once the boundary conditions are applied, the first and last element of the force vector go to 0

Non-symbolic representations of the mass and stiffness matrices and force vector with the above user input

```
[[0.03333333 0.01666667 0 0 0 0 0 0 0 0 0]
 [0.01666667 0.06666667 0.01666667 0 0 0 0 0 0 0 0]
 [0 0.01666667 0.06666667 0.01666667 0 0 0 0 0 0 0]
 [0 0 0.01666667 0.06666667 0.01666667 0 0 0 0 0 0]
 [0 0 0 0.01666667 0.06666667 0.01666667 0 0 0 0 0]
 [0 0 0 0 0.01666667 0.06666667 0.01666667 0 0 0 0]
 [0 0 0 0 0 0.01666667 0.06666667 0.01666667 0 0 0]
 [0 0 0 0 0 0 0.01666667 0.06666667 0.01666667 0 0]
 [0 0 0 0 0 0 0 0.01666667 0.06666667 0.01666667 0]
 [0 0 0 0 0 0 0 0 0.01666667 0.06666667 0.01666667]
 [0 0 0 0 0 0 0 0 0 0.01666667 0.03333333]]

 [[ 10. -10.   0.   0.   0.   0.   0.   0.   0.   0.   0.]
 [-10.  20. -10.   0.   0.   0.   0.   0.   0.   0.   0.]
 [  0. -10.  20. -10.   0.   0.   0.   0.   0.   0.   0.]
 [  0.   0. -10.  20. -10.   0.   0.   0.   0.   0.   0.]
 [  0.   0.   0. -10.  20. -10.   0.   0.   0.   0.   0.]
 [  0.   0.   0.   0. -10.  20. -10.   0.   0.   0.   0.]
 [  0.   0.   0.   0.   0. -10.  20. -10.   0.   0.   0.]
 [  0.   0.   0.   0.   0.   0. -10.  20. -10.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0. -10.  20. -10.   0.]
 [  0.   0.   0.   0.   0.   0.   0.   0. -10.  20. -10.]
 [  0.   0.   0.   0.   0.   0.   0.   0.   0. -10.  10.]]

 [[0], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.1*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.2*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.3*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.4*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.6*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.7*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.8*pi)], 
  [0.1*(-1 + pi**2)*exp(-t)*sin(0.9*pi)], 
  [0]]
 ```

## Finding a Numerical Solution
After going through the desired number of timesteps the final u_n vector is printed out for the user to see:

### Forward Euler Method

```
[0 0.0750985665896654 0.176930981481778 0.257851201568156
 0.309845339376259 0.327767761614359 0.309845339376260 0.257851201568156
 0.176930981481779 0.0750985665896655 1.22464679914735e-16]
```
![alt text](https://github.com/Phart1226/COE352/blob/master/Project2/FEGraph.JPG?raw=true)

### Backward Euler Method
```
[0 0.114716515595626 0.218203779367789 0.300331736900306 0.353061131490766 0.371230442609360 0.353061131490766 0.300331736900306 0.218203779367789 0.114716515595626 1.22464679914735e-16]
```
![alt text](https://github.com/Phart1226/COE352/blob/master/Project2/BEGraph.JPG?raw=true)

A graph for the analytical solution vs the numerical solution is printed 
for the user



## Project Questions

1. **At what dt does instability occur and how does the solution change as the number of nodes decreases?**

The 1D Galerkin method used, is unstable between dt's of 1/551 and 1/600. As the number of nodes decreases, the solution becomes less accurate and does not align with the analytical solution very well. 

2. **What happens as the timestep becomes greater than or equal to the spatial step size?** 

As the spatial step size becomes greater than or equal to the timestep, more information is forward propagating spatially than information dependent on time is, so the model is going to be inaccuarate since spatial accuracy is lost.






