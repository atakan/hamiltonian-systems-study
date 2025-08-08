import numpy as np
import matplotlib.pyplot as plt

'''"simple" integration routine for Henon-Heiles problem.
'''

# I guess I should make the following part of a class to make it more tidy, but postponing it. --ato 2025-08-07 15:54
cb2 = np.cbrt(2.0)
w0 = -cb2/(2-cb2)
w1 = 1/(2-cb2)
TJc1 = TJc4 = w1/2
TJc2 = TJc3 = (w0+w1)/2
TJd1 = TJd3 = w1
TJd2 = w0


def drift(q, p, dt):
    '''change the positions, without changing the momenta'''
    qn = np.asarray([q[0] + dt*p[0],
                     q[1] + dt*p[1]])
    pn = p.copy()
    return qn, pn

def kick(q, p, dt):
    '''change the momenta, without changing the positions'''
    qn = q.copy()
    pn = np.asarray([ p[0] + dt*(-q[0] - 2*q[0]*q[1]),
                      p[1] + dt*(-q[1] - (q[0]*q[0] - q[1]*q[1])) ])
    return qn, pn

def DKD_leapfrog(q, p, dt):
    q,p = drift(q, p, dt/2)
    q,p = kick(q, p, dt/2)
    q,p = drift(q, p, dt/2)
    return q,p

def triple_jump(q, p, dt):
    q, p = drift(q, p, dt*TJc1)
    q, p = kick(q, p, dt*TJd1)
    q, p = drift(q, p, dt*TJc2)
    q, p = kick(q, p, dt*TJd2)
    q, p = drift(q, p, dt*TJc3)
    q, p = kick(q, p, dt*TJd3)
    q, p = drift(q, p, dt*TJc4)
    return q, p

def init_cond_dankowicz(H):
    '''we plot the intersection with x=0 plane'''
    x = 0.0

    '''q1 and p1 lie in a circle of radius sqrt(2H)
       here, we choose them somewhat randomly.'''
    y = 0.3
    vy = -0.14

    vx = np.sqrt(2*H - (x*x + y*y -2*x*x*y + 2/3*y*y*y - vy*vy))

    print("initial conditions:")
    print(np.asarray([x, y]), np.asarray([vx, vy]))

    return np.asarray([x, y]), np.asarray([vx, vy])

def calc_E(q, p):
    return (0.5*(q[0]*q[0] + q[1]*q[1]) +
            (q[0]*q[0]*q[1] - q[1]*q[1]*q[1]/3) +
            (p[0]*p[0]+p[1]*p[1]))

def main():
    tend = 1240.0
    t = 0.0
    dt = 0.01

    H0 = 0.083333
    
    q_prev, p_prev = init_cond_dankowicz(H0)

    xplt = []
    yplt = []
    while t<tend:
        q, p = triple_jump(q_prev, p_prev, dt)
        #print("%8.3e %8.3e %8.3e %8.3e" % (q[0], p[0], q[1], p[1]))
        #print("     %8.3e %8.3e %8.3e" % (q_prev[1], q[1], q_prev[1]*q[1]))
        if q_prev[0]*q[0] <= 0 and p[0]>0: # x passed through zero
            xplt.append((q_prev[1]+q[1])/2) # we can be more fancy here and make linear interpolation, let's postpone though
            yplt.append((p_prev[1]+p[1])/2)
            #print("********************", xplt[-1], yplt[-1])
        q_prev = q.copy()
        p_prev = p.copy()
        t += dt
        print("%.4e" % (calc_E(q, p)-H0))
    
    # Example Matplotlib plot
    plt.figure(figsize=(8, 6))
    plt.plot(xplt, yplt, ".")
    plt.title('Henon and Heiles  figure 3')
    plt.xlabel('q1')
    plt.ylabel('p1')
    plt.grid(False)
    #plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
