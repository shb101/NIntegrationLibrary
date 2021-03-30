import ctypes
import math
import time

def NIntTrapz(func, x0, x1, N):
    engine = ctypes.cdll.LoadLibrary('NIntegrationEngine.dll')
    dx = (x1-x0)/N
    F = [func(x0+i*dx) for i in range(N+1)]
    ctypes_F = (ctypes.c_double * (N+1))(*F)
    engine.NIntTrapz.restype = ctypes.c_double
    I = engine.NIntTrapz(ctypes.c_double(x0), ctypes.c_double(x1), ctypes_F, ctypes.c_int(N))
    return I

def PyNIntTrapz(func, x0, x1, N):
    dx = (x1 - x0)/N
    I = 0
    for i in range(N+1):
        I += dx * (func(x0+dx*i) + func(x0+dx*(i+1)))/2
    return I

# Testing
#func = lambda x: math.sin(math.pi*x/2)
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