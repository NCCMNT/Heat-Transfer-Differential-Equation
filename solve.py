from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import numpy as np

def main(n):

    def k(x):
        return 1 if x <= 1 else 2 * x
    
    # create evenly spaced samples between boundries of omega
    # omega = < 0, 2 >
    samples = np.linspace(0,2,n+1)
    h = samples[1] - samples[0]

    # initialize matrices
    # B - bilinear matrix n x n
    # L - linear matrix 1 x n
    B = lil_matrix((n+1, n+1))
    L = np.zeros(n+1)

    for i in range(n):
        xi, xi1 = samples[i], samples[i+1]
        ke = k( (xi + xi1) / 2) / h
        B[i,i] += ke
        B[i, i+1] -= ke
        B[i+1, i] -= ke
        B[i+1, i+1] += ke

        fe = h * (100 * (xi + xi1) / 2) / 2
        L[i] += fe
        L[i+1] += fe

    # 1st boundary conditions
    B[0,0] -= 1  # corresponds to -v(0)u(0) in B
    L[0] -= 20   # corresponds to -20v(0) in L
    
    # 2nd boundaty condition: u(2) = 0
    B[-1,-1] = 1
    L[-1] = 0
        
    # Solving matrix equation
    B = B.tocsr()  # convert matrix to Compressed Sparse Row (CSR) for more effective computing
    u = spsolve(B, L)

    # Results
    print("Solved u:", u)

    # Visualization
    import matplotlib.pyplot as plt
    plt.plot(samples, u, label='u(x)')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Heat transfer equation')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    while True:
        n = int(input("Set number of samples: "))
        if n >= 2: break
        print("Incorrect value")
    main(n)