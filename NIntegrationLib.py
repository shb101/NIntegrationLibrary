import ctypes
import math
import time

engine = ctypes.cdll.LoadLibrary('NIntegrationEngine.dll')

def NIntTrapz(func, x0, x1, N):
    dx = (x1-x0)/N
    F = [func(x0+i*dx) for i in range(N+1)]
    c_F = (ctypes.c_double * (N+1))(*F)
    engine.NIntTrapz.restype = ctypes.c_double
    I = engine.NIntTrapz(ctypes.c_double(x0), ctypes.c_double(x1), c_F, ctypes.c_int(N))
    return I

def NIntTrapz2D(func, x0, x1, y0, y1, N, M):
    dx = (x1 - x0)/N
    c_F = (ctypes.POINTER(ctypes.c_double) * (N+1))()
    c_y_lb = (ctypes.c_double * (N+1))()
    c_y_ub = (ctypes.c_double * (N+1))()
    for i in range(N+1):
        x = x0 + i*dx
        lb = y0(x)
        ub = y1(x)
        c_y_lb[i] = lb
        c_y_ub[i] = ub
        dy = (ub - lb)/M
        c_F[i] = (ctypes.c_double * (M+1))(*[func(x, lb+dy*j) for j in range(M+1)])
    engine.NIntTrapz2D.restype = ctypes.c_double
    I = engine.NIntTrapz2D(ctypes.c_double(x0), ctypes.c_double(x1), c_y_lb, c_y_ub, c_F, ctypes.c_int(N), ctypes.c_int(M))
    return I


def PyNIntTrapz(func, x0, x1, N):
    dx = (x1 - x0)/N
    I = 0
    for i in range(N+1):
        I += dx * (func(x0+dx*i) + func(x0+dx*(i+1)))/2
    return I

# Testing
def test_func1():
    func = lambda x: math.exp(-x**2/2) * math.sin(math.sqrt(x))
    N = 1000000

    print('1. C++')
    t1 = time.time()
    print(NIntTrapz(func,0,1,N))
    t2 = time.time()
    print('Time taken = ' + str(t2-t1) + ' seconds')

    print('2. Py')
    t1 = time.time()
    print(PyNIntTrapz(func,0,1,N))
    t2 = time.time()
    print('Time taken = ' + str(t2-t1) + ' seconds')

def main():
    func = lambda x, y: math.cos(x**2 * y) * math.exp(-(x**2 + y**2)/4)
    x0 = 0
    x1 = 1
    y0 = lambda x: x
    y1 = lambda x: 2*x
    N = 1000
    M = 1000
    I = NIntTrapz2D(func, x0, x1, y0, y1, N, M)
    print(I)

if __name__ == '__main__':
    main()
