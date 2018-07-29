import numpy as np

def change(*argv):
    for arg in argv:
        arg[3] = 3


a = np.random.randn(10)
print a
change(a[3:])
print a
