from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

def main(n):
    h = 2 / n
    # create evenly spaced samples between boundries of omega
    # omega = < 0, 2 >
    samples = np.linspace(0, 2, n + 1)
    # h = samples[1] - samples[0]

    def k(x):
        if x < 0 or x > 2:
            return 1
        if x <= 1:
            return 1
        return 2 * x
    
    def e(i, x):
        xi = i * h
        xi1 = (i + 1) * h
        if xi <= x < xi1:
            return (x - xi) / h
        xi_1 = (i - 1) * h
        if xi_1 <= x < xi:
            return (xi - x) / h
        return 0

    def de_dx(i, x):
        xi = i * h
        xi1 = (i + 1) * h
        if xi <= x < xi1:
            return 1 / h
        xi_1 = (i - 1) * h
        if xi_1 <= x < xi:
            return -1 / h
        return 0

    # Bilinear form B
    def B(i, j):
        if abs(i - j) > 1:
            return 0
        a = max(0, (i - 1) * h) # x_(i-1)
        b = min(2, (i + 1) * h) # x_(i+1)
        integral = integrate.quad(lambda x: k(x) * de_dx(i, x) * de_dx(j, x), a, b)[0] - e(i, 0) * e(j, 0)
        return integral
    
    # Linear form L
    def L(i):
        a = max(0, (i - 1) * h) # x_(i-1)
        b = min(2, (i + 1) * h) # x_(i+1)
        return (integrate.quad(lambda x: 100 * x * e(i, x), a, b)[0] - 20 * e(i, 0))
    
    def show_plot(x,y):
        # Visualization
        plt.plot(x, y, label='u(x)')
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.title(f'Heat transfer equation - {n} samples')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    # initialize matrices
    # B - bilinear matrix n x n
    # L - linear matrix 1 x n
    B_matrix = lil_matrix((n + 1, n + 1))
    L_vector = np.zeros(n + 1)

    # fill matrices with values
    for i in range(n + 1):
        for j in range(max(0, i - 1), min(n + 1, i + 2)):
            B_matrix[i, j] = B(i, j)
        L_vector[i] = L(i)
    
    # 1st boundary conditions
    B_matrix[0, :] = 0  
    B_matrix[0, 0] = 1
    L_vector[0] = -20 # corresponds to -20v(0) in L

    # 2nd boundaty condition: u(2) = 0
    B_matrix[-1, :] = 0
    B_matrix[-1, -1] = 1
    L_vector[-1] = 0

    # Solve matrix equation
    B_matrix = B_matrix.tocsr() # convert matrix to Compressed Sparse Row (CSR) for more effective computing
    u = spsolve(B_matrix, L_vector)
    
    show_plot(samples, u)

if __name__ == '__main__':
    while True:
        n = int(input("Set number of samples: "))
        if n >= 2: break
        print("Incorrect value")
    main(n)