from subprocess import Popen, PIPE
import math

def NIntTrapz(func, x0, x1, N):
    p = Popen(['NIntEngine.exe'], shell=True, stdout=PIPE, stdin=PIPE)
    # Input N, x0, x1:
    in1 = bytes('0 ' + str(N) + ' ' + str(x0) + ' ' + str(x1) + '\n', 'UTF-8')
    p.stdin.write(in1)
    p.stdin.flush()
    # Input func values:
    dx = (x1-x0)/N
    for i in range(N+1):
        in2 = bytes(str(func(x0 + i*dx)) + '\n', 'UTF-8')
        p.stdin.write(in2)
        p.stdin.flush()
    # Read the result:
    result = p.stdout.readline().strip()
    return result

# Testing
func = lambda x: math.sin(math.pi*x/2)
print(NIntTrapz(func,0,1,10000))