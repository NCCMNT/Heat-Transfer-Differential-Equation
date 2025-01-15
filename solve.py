from scipy import integrate

def solve(n):
    h = 1 / n

    def e(i,x):
        x_im = (i - 1) * h
        x_i = i * h
        x_ip = (i + 1) * h

        if x_im < x < x_i:
            return (x - x_im) / (x_i - x_im)
        if x_i < x < x_ip:
            return (x_ip - x) / (x_ip - x_i)
        return 0
    
    def ediff(i,x):
        x_im = (i - 1) * h
        x_i = i * h
        x_ip = (i + 1) * h

        if x_im < x < x_i:
            return n
        if x_i < x < x_ip:
            return -n
        return 0

    def k(x):
        if 0 <= x <= 1:
            return 1
        if 1 < x <= 2:
            return 2*x
        return 0
    
    def B(i, j, a, b):
        return -(e(i,0) * e(j,0)) * integrate.quad(lambda x: k(x) * ediff(i,x) * ediff(j,x), a, b)[0]
    
    def L(i, a, b):
        return integrate.quad(lambda x: 100 * x * e(i,x), a, b) - 20 * e(i,0)
    
    def matrixB():
        pass
