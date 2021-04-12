import ctypes
import math
import time

# Load the numerical integration engine:
ENGINEPATH = 'NIntegrationEngine/NIntegrationEngine.dll'
engine = ctypes.cdll.LoadLibrary(ENGINEPATH)
# Settings
engine.NIntTrapz.restype = ctypes.c_double
engine.NIntTrapz2D.restype = ctypes.c_double
engine.NIntTrapz3D.restype = ctypes.c_double
engine.NIntSimpson.restype = ctypes.c_double
engine.NIntSimpson2D.restype = ctypes.c_double
engine.NIntSimpson3D.restype = ctypes.c_double

# 1D Numerical Integration
def NInt(func, x0, x1, N, method='trapz'):
    dx = (x1-x0)/N
    c_F = (ctypes.c_double * (N+1))()
    for i in range(N+1):
        c_F[i] = func(x0 + i*dx)
    if method == 'trapz':
        I = engine.NIntTrapz(ctypes.c_double(x0), ctypes.c_double(x1), c_F, ctypes.c_int(N))
    elif method == 'simpson':
        I = engine.NIntSimpson(ctypes.c_double(x0), ctypes.c_double(x1), c_F, ctypes.c_int(N))
    return I

# 2D Numerical Integration
def NInt2D(func, x0, x1, y0, y1, N, M, method='trapz'):
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
        c_F[i] = (ctypes.c_double * (M+1))()
        for j in range(M+1):
            c_F[i][j] = func(x, lb + j*dy)
    if method == 'trapz':
        I = engine.NIntTrapz2D(ctypes.c_double(x0), ctypes.c_double(x1), c_y_lb, c_y_ub, c_F, ctypes.c_int(N), ctypes.c_int(M))
    elif method == 'simpson':
        I = engine.NIntSimpson2D(ctypes.c_double(x0), ctypes.c_double(x1), c_y_lb, c_y_ub, c_F, ctypes.c_int(N), ctypes.c_int(M))
    return I

# 3D Numerical Integration
def NInt3D(func, x0, x1, y0, y1, z0, z1, N, M, L, method='trapz'):
    dx = (x1 - x0)/N
    c_F = (ctypes.POINTER(ctypes.POINTER(ctypes.c_double)) * (N+1))()
    c_y_lb = (ctypes.c_double * (N+1))()
    c_y_ub = (ctypes.c_double * (N+1))()
    c_z_lb = (ctypes.POINTER(ctypes.c_double) * (N+1))()
    c_z_ub = (ctypes.POINTER(ctypes.c_double) * (N+1))()
    for i in range(N+1):
        x = x0 + i*dx
        c_y_lb[i] = y0(x)
        c_y_ub[i] = y1(x)
        dy = (c_y_ub[i] - c_y_lb[i])/M
        c_F[i] = (ctypes.POINTER(ctypes.c_double) * (M+1))()
        c_z_lb[i] = (ctypes.c_double * (M+1))()
        c_z_ub[i] = (ctypes.c_double * (M+1))()
        for j in range(M+1):
            y = c_y_lb[i] + j*dy
            c_z_lb[i][j] = z0(x,y)
            c_z_ub[i][j] = z1(x,y)
            dz = (c_z_ub[i][j] - c_z_lb[i][j])/L
            c_F[i][j] = (ctypes.c_double * (L+1))()
            for k in range(L+1):
                c_F[i][j][k] = func(x,y,c_z_lb[i][j]+k*dz)
    if method == 'trapz':
        I = engine.NIntTrapz3D(ctypes.c_double(x0), ctypes.c_double(x1), c_y_lb, c_y_ub, c_z_lb, c_z_ub, c_F, ctypes.c_int(N), ctypes.c_int(M), ctypes.c_int(L))
    elif method == 'simpson':
        I = engine.NIntSimpson3D(ctypes.c_double(x0), ctypes.c_double(x1), c_y_lb, c_y_ub, c_z_lb, c_z_ub, c_F, ctypes.c_int(N), ctypes.c_int(M), ctypes.c_int(L))
    return I

def PyNIntTrapz(func, x0, x1, N):
    dx = (x1 - x0)/N
    I = 0
    F = [func(x0 + i*dx) for i in range(N+1)]
    for i in range(N):
        I += dx * (F[i] + F[i+1])/2
    return I

def PyNIntTrapz2D(func, x0, x1, y0, y1, N, M):
    dx = (x1 - x0)/N
    I = 0
    F = []
    y_0 = []
    y_1 = []
    dy = []
    for i in range(N+1):
        F.append([])
        x = x0 + i*dx
        y_0.append(y0(x))
        y_1.append(y1(x))
        dy.append((y_1[i] - y_0[i])/M)
        for j in range(M+1):
            F[i].append(func(x, y_0[i] + j*dy[i]))
    for i in range(N):
        for j in range(M):
            I += dx*dy[i]*(F[i][j] + F[i+1][j] + F[i][j+1] + F[i+1][j+1])/4
    return I

# Testing
def test_func1():
    func = lambda x: math.exp(-x**2/2) * math.sin(math.sqrt(x))
    N = 1000000

    print('1) C')
    t1 = time.time()
    print(NInt(func,0,1,N, method='simpson'))
    t2 = time.time()
    print('Time taken = ' + str(t2-t1) + ' seconds')

    print('2) Py')
    t1 = time.time()
    print(PyNIntTrapz(func,0,1,N))
    t2 = time.time()
    print('Time taken = ' + str(t2-t1) + ' seconds')

def test_func2():
    func = lambda x, y: math.cos(x**2 * y) * math.exp(-(x**2 + y**2)/4)
    x0 = 0
    x1 = 1
    y0 = lambda x: x
    y1 = lambda x: 2*x
    N = 4000
    M = 4000
    # C
    print('1) C Code')
    t1 = time.time()
    I = NInt2D(func, x0, x1, y0, y1, N, M, method='simpson')
    t2 = time.time()
    print(I)
    print('Time taken = ' + str(t2-t1) + ' seconds')
    # Py
    print('2) Py Code')
    t1 = time.time()
    I = PyNIntTrapz2D(func, x0, x1, y0, y1, N, M)
    t2 = time.time()
    print(I)
    print('Time taken = ' + str(t2-t1) + ' seconds')

def main():
    func = lambda x, y, z: math.sin(x**2 * y * z) * math.sqrt(x**2 + y**2 + z**2)
    x0 = 0
    x1 = 1
    y0 = lambda x: 0.1*x
    y1 = lambda x: 0.7*x
    z0 = lambda x, y: x - y
    z1 = lambda x, y: x + y
    N = 100
    M = 100
    L = 100
    print(NInt3D(func, x0, x1, y0, y1, z0, z1, N, M, L))

if __name__ == '__main__':
    test_func2()
