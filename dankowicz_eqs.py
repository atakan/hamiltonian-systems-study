import numpy as np

def dqdt(q, p, t):
    return np.asarray([
        p[0],
        p[1],
    ])

def dpdt(q, p, t):
    return np.asarray([
        -2*q[0]*q[1] - q[0],
        -q[0]**2 - 2*q[1],
    ])

