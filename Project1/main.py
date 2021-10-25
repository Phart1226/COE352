import numpy as np

def solve_disp(K, f):
    try:
        u = np.linalg.solve(K, f*9.81)
    except np.linalg.LinAlgError as err:
        if 'Singular matrix' in str(err):
            print('The stiffness matrix K is singluar and therefore not invertible, a psuedo inverse calculation to give a close approximation to a real solution will be calculated\n')
            u = np.matmul(np.linalg.pinv(K), f*9.81)
    return u

def solve_stress(C, e):
    return np.matmul(C, e)

def solve_elong(A, u):
    return np.matmul(A, u)

def calc_svd(A):
    U, s, VH = np.linalg.svd(A)
    return s

def find_sing_eig(A, C, A_t, K):
    sing_A = calc_svd(A)
    sing_C = calc_svd(C)
    sing_A_t = calc_svd(A_t)
    sing_K = calc_svd(K)

    print('A matrix: ')
    print(A)
    print('Singular values of A')
    print(sing_A)
    print('Eigenvalues of A')
    print(np.square(sing_A))
    print()
    print('Condition number of A')
    print("{:.3f}".format(calc_cond(sing_A)))
    print()
    
    print('C matrix: ')
    print(C)
    print('Singular values of C')
    print(sing_C)
    print('Eigenvalues of C')
    print(np.square(sing_C))
    print()
    print('Condition number of C')
    print("{:.3f}".format(calc_cond(sing_C)))
    print()

    print('A_t matrix: ')
    print(A_t)
    print('Singular values of A_t')
    print(sing_A_t)
    print('Eigenvalues of A_t')
    print(np.square(sing_A_t))
    print()
    print('Condition number of A_t')
    print("{:.3f}".format(calc_cond(sing_A_t)))
    print()

    print('Stiffness matrix K: ')
    print(K)
    print('Singular values of K')
    print(sing_K)
    print('Eigenvalues of K')
    print(np.square(sing_K))
    print()
    print('Condition number of K')
    print("{:.3f}".format(calc_cond(sing_K)))
    print()

def calc_cond(sing_vals):
    return abs(sing_vals[0])/abs(sing_vals[-1])

def create_A_t(condition, num_masses, num_springs):
    A_t = np.zeros((num_masses, num_springs))

    # fixed-free scenario
    if condition == 0:
        spot = 0
        for row in range(num_masses):
            A_t[row,spot] = 1
            spot += 1
            if spot != num_masses:
                A_t[row,spot] = -1

    # fixed-fixed scenario
    elif condition == 1:
        A_t = np.zeros((num_masses, num_springs))
        spot = 0
        for row in range(num_masses):
            A_t[row,spot] = 1
            spot += 1
            A_t[row,spot] = -1
    
    # free-free scenario
    else:
        A_t = np.zeros((num_springs, num_masses))
        spot = 0
        for row in range(num_springs):
            A_t[row,spot] = -1
            spot += 1
            A_t[row,spot] = 1
        A_t = A_t.transpose()

    return A_t

def create_C (spring_const):
    return  np.diag(spring_const)

def generate_k_mat(spring_const, condition, num_masses):
    num_springs = len(spring_const)
    A_t = create_A_t(condition, num_masses, num_springs)  
    C = create_C(spring_const)
    A = A_t.transpose()
    K = np.zeros((num_masses, num_springs))

    

    if condition == 0:
        K = np.matmul(np.matmul(A_t,C), A)
    else:
        K = np.matmul(A_t, np.matmul(C, A))

    find_sing_eig(A,C,A_t,K)

    return K



def main ():
    # gather user input
    print('Welcome to the Mass-Spring Simulator')
    print("Choose either Fixed-Free, Fixed-Fixed, or Free-Free boundary condition (enter the number of the condition): ")
    print("[0]: Fixed Free")
    print("[1]: Fixed Fixed")
    print("[-1]: Free Free")
    print()
    condition = int(input("Boundary condition: "))
    num_masses = int(input('Enter the number of masses: '))
    masses = np.zeros(num_masses)
    for mass in range(num_masses):
        print('Enter weight of mass {}:'.format(mass + 1))
        masses[mass] = int(input())
    print()

    spring_const = np.zeros(num_masses + condition)
    for spring in range(len(spring_const)):
        print('Enter spring constant {}:'.format(spring + 1))
        spring_const[spring] = int(input())
    print()

    # generate k matrix
    K = generate_k_mat(spring_const, condition, len(masses))
    
    
    # solve for displacement
    u = solve_disp(K, masses)
    for disp in range(len(u)):
        print('The displacement of mass ' + str(disp + 1) + ' is ' + "{:.2f}".format(u[disp]))
    print()

    # solve for elongation
    e = solve_elong(create_A_t(condition, len(masses), len(spring_const)).transpose(), u)
    for elong in range(len(e)):
        print('The elongation of spring ' + str(elong + 1) + ' is ' + "{:.2f}".format(e[elong]))
    print()

    # solve for internal stress
    w = solve_stress(create_C(spring_const), e)
    for stress in range(len(w)):
        print('The internal stress of spring ' + str(stress + 1) + ' is ' + "{:.2f}".format(w[stress]))
    print()

    return




if __name__ == "__main__":
    main()