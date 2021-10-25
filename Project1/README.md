---
# Project 1  Mass-Spring System Simulator

### Preston Hart


#### Objective

To create a software that calculates the displacement of the masses and the elongation and internal stresses of the springs in a hanging mass spring system at its equilibrium state for different boundary conditions.

#### How to Use

1. After downloading the zip folder, make sure that numpy has been pip installed in the folder that the software will run in.

2. Navigate to the main.py file and run it in a command window using the command:

```
python3 main.py
```

3. The mass-spring simulator software will appear and the user will enter in the their boundary condition and number of masses

#### Example input/output

```
\COE 352>python main.py
Welcome to the Mass-Spring Simulator
Choose either Fixed-Free, Fixed-Fixed, or Free-Free boundary condition (enter the number of the condition):
[0]: Fixed Free
[1]: Fixed Fixed
[-1]: Free Free

Boundary condition: 0
Enter the number of masses: 3
Enter weight of mass 1:
1
Enter weight of mass 2:
1
Enter weight of mass 3:
1

Enter spring constant 1:
1
Enter spring constant 2:
1
Enter spring constant 3:
1

A matrix:
[[ 1.  0.  0.]
 [-1.  1.  0.]
 [ 0. -1.  1.]]
Singular values of A
[1.80193774 1.2469796  0.44504187]
Eigenvalues of A
[3.2469796  1.55495813 0.19806226]

Condition number of A
4.049

C matrix:
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
Singular values of C
[1. 1. 1.]
Eigenvalues of C
[1. 1. 1.]

Condition number of C
1.000

A_t matrix:
[[ 1. -1.  0.]
 [ 0.  1. -1.]
 [ 0.  0.  1.]]
Singular values of A_t
[1.80193774 1.2469796  0.44504187]
Eigenvalues of A_t
[3.2469796  1.55495813 0.19806226]

Condition number of A_t
4.049

Stiffness matrix K:
[[ 2. -1.  0.]
 [-1.  2. -1.]
 [ 0. -1.  1.]]
Singular values of K
[3.2469796  1.55495813 0.19806226]
Eigenvalues of K
[10.54287655  2.41789479  0.03922866]

Condition number of K
16.394

The displacement of mass 1 is 29.43
The displacement of mass 2 is 49.05
The displacement of mass 3 is 58.86

The elongation of spring 1 is 29.43
The elongation of spring 2 is 19.62
The elongation of spring 3 is 9.81

The internal stress of spring 1 is 29.43
The internal stress of spring 2 is 19.62
The internal stress of spring 3 is 9.81
```

#### The Free-Free Boundary Condition

In reality a mass-spring system cannot have a free-free boundary condition, since masses and springs cannot be connected and yet just float in space. Therefore, a solution to a system with this specific boundary condition has no real-world implications or meaning.

However, if the general mathmatical procudure for finding a solution to the mass-spring system is followed for this boundary condition, then a solution is still obtainable. 

From the elongation equations, the system is undertermined for two free ends becuase there is one more mass than springs and no boundary conditions to apply. This poses a problem when solving for the displacements becuase the stiffness matrix, K, is not invertible. However, the psuedo-inverse of the K matrix can be used to solve the system.

Because the masses can move without stretching or compressing the springs in a free-free system, it imposes difficulty for the math, since the displacements of the masses are related to the springs' elongations. So this can explain why the outputs for a system with two free ends do not always make since. 