import numpy as np

def dqdt(q, p, t):
    return np.asarray([
        p[0],
        p[1],
    ])

def dpdt(q, p, t):
    return np.asarray([
        - q[0],
        - 2*q[1],
    ])

