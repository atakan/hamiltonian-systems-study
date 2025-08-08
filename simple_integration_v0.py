import numpy as np
import matplotlib.pyplot as plt

#from dankowicz_eqs import dqdt, dpdt
from dankowicz_eqs_simp import dqdt, dpdt

def drift(q, p, dt):
    '''change the positions, without changing the momenta'''
    return q+dqdt(q, p, 0)*dt, p

def kick(q, p, dt):
    '''change the momenta, without changing the positions'''
    return q, p + dpdt(q, p, 0)*dt

def DKD_leapfrog(q, p, dt):
    q,p = drift(q, p, dt/2)
    q,p = kick(q, p, dt/2)
    q,p = drift(q, p, dt/2)
    return q,p

def init_cond_dankowicz(H):
    '''we plot the intersection with q2=0 plane'''
    q2 = 0.0

    '''q1 and p1 lie in a circle of radius sqrt(2H)
       here, we choose them somewhat randomly.'''
    q1 = 0.4*np.sqrt(2*H)
    p1 = -0.15*np.sqrt(2*H)

    p2 = np.sqrt(2*H - (p1*p1 + q1*q1 + 2*q2*q2 - 2*q1*q1*q2))

    print("initial conditions:")
    print(np.asarray([q1, q2]), np.asarray([p1, p2]))

    return np.asarray([q1, q2]), np.asarray([p1, p2])

def main():
    tend = 124.0
    t = 0.0
    dt = 0.01

    H = 0.62
    
    q_prev, p_prev = init_cond_dankowicz(H)

    xplt = []
    yplt = []
    while t<tend:
        q, p = DKD_leapfrog(q_prev, p_prev, dt)
        print("%8.3e %8.3e %8.3e %8.3e" % (q[0], p[0], q[1], p[1]))
        print("     %8.3e %8.3e %8.3e" % (q_prev[1], q[1], q_prev[1]*q[1]))
        if q_prev[1]*q[1] <= 0: # q2 passed through zero
            xplt.append((q_prev[0]+q[0])/2) # we can be more fancy here and make linear interpolation, let's postpone though
            yplt.append((p_prev[0]+p[0])/2)
            print("********************", xplt[-1], yplt[-1])
        q_prev = q.copy()
        p_prev = p.copy()
        t += dt
    
    # Example Matplotlib plot
    plt.figure(figsize=(8, 6))
    plt.plot(xplt, yplt, ".")
    plt.title('Dankowicz figure on p35')
    plt.xlabel('q1')
    plt.ylabel('p1')
    plt.grid(False)
    #plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
