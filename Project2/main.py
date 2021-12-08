import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def eval_force_vec(F, ctime, nodes_arr):
    F = F.subs('t', ctime)
    return [F[n].subs('X', nodes_arr[n]).evalf() for n in range(len(nodes_arr))]

def euler_calc(method, u_n, nodes_arr, dt, M, K, F, bc_arr):
    B = (1/dt)*M + K
    I = np.eye(len(u_n))
    for n in range(int(1/dt)):
        ctime = n*dt  # current time
        # Forward Euler
        if method == 'FE':
            iden = (I - dt*(np.linalg.inv(M)@K))
            # apply boundary conditions of (0,1)x(0,1)
            iden[0] = bc_arr[0]
            iden[-1] = bc_arr[1]
            iden [0][0] = bc_arr[2]
            iden[-1][-1] = bc_arr[3]
            update_force = dt*(np.linalg.inv(M)@eval_force_vec(F, ctime, nodes_arr))
            # apply boundary conditions of (0,1)x(0,1)
            update_force[0] = bc_arr[4]
            update_force[-1] = bc_arr[5]

            u_n = iden@u_n + np.array(update_force)
        # Backward Euler
        else:
            part1 = (1/dt)*(np.linalg.inv(B)@M)
            # apply boundary conditions of (0,1)x(0,1)
            part1[0] = bc_arr[0] 
            part1[-1] = bc_arr[1]
            part1 [0][0] = bc_arr[2]
            part1[-1][-1] = bc_arr[3]
            part2 = np.linalg.inv(B)@eval_force_vec(F, ctime, nodes_arr)
            # apply boundary conditions of (0,1)x(0,1)
            part2[0] = bc_arr[4]
            part2[-1] = bc_arr[5]

            u_n = part1@u_n + part2
            
    return u_n
def phi_r(r, X, deg):
    if isinstance(X, sp.Symbol):
        h = sp.Rational(1, deg)  # node spacing
        nodes = [2*i*h - 1 for i in range(deg+1)]
    else:
        nodes = np.linspace(-1, 1, deg+1)
    return Lagrange_polynomial(X, r, nodes)

def Lagrange_polynomial(x, i, points):
    p = 1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p

def basis(deg=1):
    X = sp.Symbol('X')
    phi = [phi_r(r, X, deg) for r in range(deg+1)]
    return phi

def element_matrix(phi, elem_loc, symbolic=True):
    n = len(phi)
    M = sp.zeros(n, n)
    K = sp.zeros(n, n)
    X = sp.Symbol('X')

    if symbolic:
        h = sp.Symbol('h')
    else:
        h = elem_loc[1] - elem_loc[0]
    detJ = h/2  # for 1D FE
    for r in range(n):
        for s in range(r, n):
            M[r,s] = sp.integrate(phi[r]*phi[s]*detJ, (X, -1, 1))
            M[s,r] = M[r,s]
            K[r,s] = sp.integrate(sp.diff(phi[r],X)*sp.diff(phi[s],X)*(1/detJ), (X, -1, 1))
            K[s,r] = K[r,s]

    return M, K

def quadrature_points(f):
    quad_points = [1/np.sqrt(3), -1/np.sqrt(3)]
    return f.subs(sp.Symbol('x'), quad_points[0]) +  f.subs(sp.Symbol('x'), quad_points[1])

def element_vector(f, phi, elem_loc, symbolic=True):
    # Make f a function of X
    X = sp.Symbol('X')
    if symbolic:
        h = sp.Symbol('h')
    else:
        h = elem_loc[1] - elem_loc[0]
    detJ = h/2  # for 1D FE
    f_e = [quadrature_points(f*phi[r]*detJ) for r in range(2)]
    
    if symbolic:
        return f_e
    else:
        for r in range(2):
            f_e[r] = f_e[r].subs(X, elem_loc[r])
        return f_e

# Constructs the Mass, Stiffness, and Force vectors from the elemental matrices
def assemble(node_arr, elem_arr, phi, f, symbolic):
    # initalize matrices and determine if they need to be a symbolic representation
    N_n, N_e = len(node_arr), len(elem_arr)
    mlocal = np.zeros((N_n,N_n))
    klocal = np.zeros((N_n,N_n))
    F = sp.zeros(N_n,1)
    if symbolic:
        M = sp.zeros(N_n, N_n)
        K = sp.zeros(N_n, N_n)
        mlocal = sp.zeros(N_n,N_n)
        klocal = sp.zeros(N_n,N_n)
        
    else:
        M = np.zeros((N_n, N_n))
        K = np.zeros((N_n, N_n))
        mlocal = np.zeros((N_n,N_n))
        klocal = np.zeros((N_n,N_n))

    for k in range(N_e):
        #construct elemental mass and stiffness matrices, and force vector
        nodes = [node_arr[elem_arr[k][0]],node_arr[elem_arr[k][-1]]]
        M_e, K_e = element_matrix(phi, nodes, symbolic) # where mapping from x-space to xi-space occurs
        f_e = element_vector(f, phi, nodes, symbolic)

        for l in range(2):
            for m in range(2):
                mlocal[l,m] = M_e[l,m]
                klocal[l,m] = K_e[l,m]

        for l in range(2):
            global_node = elem_arr[k][l]
            for m in range(2):
                global_node2 = elem_arr[k][m]
                M[global_node,global_node2] += mlocal[l,m]
                K[global_node,global_node2] += klocal[l,m]
            F[global_node]+= f_e[l]

    return M,K,F

def plot_num_sol(u, nodes_arr, u_n, dt, method):
    plt.plot(nodes_arr, u, label='Analytical')
    plt.plot(nodes_arr, u_n, label='Numerical')
    plt.xlabel('Node location(x)')
    plt.ylabel('y')
    plt.title(method + ' Comparison of Analytical and Numerical Solutions (iterations=' + dt + ')' )
    plt.legend()
    plt.show()


def main ():
    x = sp.Symbol('X')
    t = sp.Symbol('t')
    f = (sp.pi*sp.pi - 1)*sp.E**(-t)*sp.sin(sp.pi*x)
    u = sp.E**(-1)*sp.sin(sp.pi*x)
    bc_arr = [] # array that holds the values of the boundary conditions
    # gather user input
    print("Welcome to the Finite Element Simulator")
    print("Please fill in the following user input:\n\n")
    num_nodes = int(input("Number of nodes: "))
    to_loc = int(input("Location of last node: "))
    for i in range(6):
        bc = int(input('Boundary condition ' + str(i+1) + ': '))

        bc_arr.append(bc)

    dt = 1/int(input("Number of timesteps: "))
    symbolic_rep = input("Symbolic representation for debugging?(True or False) ")
    method = input("FE or BE? ")
    if symbolic_rep == "True":
        symbolic_rep = True
    else:
        symbolic_rep = False
    
    # assemble element and element location arrays
    nodes_arr = np.linspace(0, to_loc, num_nodes).tolist() # array starting location of each node
    elem_arr = [[e + i for i in range(2)] \
                for e in range(num_nodes - 1)] # connectivity array, shows the two nodes that an element is between

    u_n = np.array([np.sin(np.pi*x_i) for x_i in nodes_arr]) # inital u_n vector creation

    M, K, F = assemble(nodes_arr, elem_arr, basis(deg=1), f, symbolic=symbolic_rep)


    u_n = euler_calc(method, u_n, nodes_arr, dt, M, K, F, bc_arr) # preforms forward or backward euler calculation for given timestep 

    u = [u.subs('X', x_i) for x_i in nodes_arr] # building analytical solution for graphing
    print(u_n) # print out of final solution vector that is graphed
    if method == 'FE':
        plot_num_sol(u, nodes_arr, u_n, str(int(1/dt)),'Forward Euler')
    else:
        plot_num_sol(u, nodes_arr, u_n, str(int(1/dt)), 'Backward Euler')
    

    return




if __name__ == "__main__":
    main()