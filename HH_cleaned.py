import numpy as np
import matplotlib.pyplot as plt

'''"simple" integration routine for Henon-Heiles problem.
this version will hopefully eat up less memory. --ato 2025-08-08 9:48
it does, but the problem was not in place vs. functional approach. So, I probably
should go back to functional approach. --ato 2025-08-09 15:27
'''

# I guess I should make the following part of a class to make it more tidy, but postponing it. --ato 2025-08-07 15:54
cb2 = np.cbrt(2.0)
TJc1 = TJc4 = 1/2/(2-cb2)
TJc2 = TJc3 = (1-cb2)*TJc1
TJd1 = TJd3 = 1/(2-cb2)
TJd2 = -cb2*TJd1

# print("# sanity check: %f %f %e %e" % (TJc1+TJc2+TJc3+TJc4,
#                                        TJc1+TJc2+TJc3+TJc4-1,
#                                        TJd1+TJd2+TJd3,
#                                        TJd1+TJd2+TJd3-1)) 

# def drift(q, p, dt):
#     '''change the positions, without changing the momenta'''
#     qn = np.asarray([q[0] + dt*p[0],
#                      q[1] + dt*p[1]])
#     pn = p.copy()
#     return qn, pn

def drift_inplace(q, p, dt):
    q[0] += dt*p[0]
    q[1] += dt*p[1]

# def kick(q, p, dt):
#     '''change the momenta, without changing the positions'''
#     qn = q.copy()
#     pn = np.asarray([ p[0] + dt*(-q[0] - 2*q[0]*q[1]),
#                       p[1] + dt*(-q[1] - (q[0]*q[0] - q[1]*q[1])) ])
#     return qn, pn

def kick_inplace(q, p, dt):
    p[0] += dt*(-q[0] - 2*q[0]*q[1])
    p[1] += dt*(-q[1] - q[0]*q[0] + q[1]*q[1])


# def DKD_leapfrog(q, p, dt):
#     q,p = drift(q, p, dt/2)
#     q,p = kick(q, p, dt/2)
#     q,p = drift(q, p, dt/2)
#     return q,p

def DKD_LF_inplace(q, p, dt):
    drift_inplace(q, p, dt/2)
    kick_inplace(q, p, dt)
    drift_inplace(q, p, dt/2)

    
# def triple_jump(q, p, dt):
#     q, p = drift(q, p, dt*TJc1)
#     q, p = kick(q, p, dt*TJd1)
#     q, p = drift(q, p, dt*TJc2)
#     q, p = kick(q, p, dt*TJd2)
#     q, p = drift(q, p, dt*TJc3)
#     q, p = kick(q, p, dt*TJd3)
#     q, p = drift(q, p, dt*TJc4)
#     return q, p


def triple_jump_inplace(q, p, dt):
    drift_inplace(q, p, dt*TJc1)
    kick_inplace(q, p, dt*TJd1)
    drift_inplace(q, p, dt*TJc2)
    kick_inplace(q, p, dt*TJd2)
    drift_inplace(q, p, dt*TJc3)
    kick_inplace(q, p, dt*TJd3)
    drift_inplace(q, p, dt*TJc4)

def init_cond_henonheiles(y, vy, H):
    '''we plot the intersection with x=0 plane'''
    x = 0.0

    vx = np.sqrt(2*H - (x*x + y*y + 2*x*x*y - 2/3*y*y*y) - vy*vy)

    print("# initial conditions:")
    print("#", np.asarray([x, y]), np.asarray([vx, vy]))

    return np.asarray([x, y]), np.asarray([vx, vy])

def calc_E(q, p):
    return (0.5*(q[0]*q[0] + q[1]*q[1] +
                 2*q[0]*q[0]*q[1] - 2/3*q[1]*q[1]*q[1]) +
            0.5*(p[0]*p[0]+p[1]*p[1]) )

def main():
    tend = 1200.0
    t = 0.0
    dt = 0.01

    H0 = 0.08333

    # y_vy = [ (0.3, -0.14), (0.31, -0.17), (0.32, -0.2)]
    # y_vy = [ (-0.2, 0.2), (-0.1, 0.2), (-0.08, 0.2) ] # basari1
    # y_vy = [ (-0.175, 0.2), (-0.15, 0.2), (-0.125, 0.2) ] # basari2
    # y_vy = [ (-0.1075, 0.2), (-0.105, 0.2), (-0.1025, 0.2) ] # basari4
    #y_vy = [ (-0.104, 0.2), (-0.1035, 0.2), (-0.103, 0.2) ] # basari5
    # y_vy = [ (-0.10475, 0.2), (-0.1045, 0.2), (-0.10425, 0.2) ] # basari6
    # y_vy = [ (-0.10495, 0.2), (-0.1049, 0.2), (-0.10485, 0.2) ] # basari7
    #y_vy = [ (-0.104825, 0.2), (-0.1048, 0.2), (-0.104775, 0.2) ] # basari8
    #y_vy = [ (-0.10481, 0.2), (-0.1048, 0.2), (-0.10479, 0.2) ] # basari9
    y_vy = [ (-0.3, 0.24),
             (-0.2, 0.2), (-0.1048, 0.2), (-0.08, 0.2),
             (-0.08, -0.2), (0.2, 0.1), (0.21113, -0.27450),
             (0.028, 0.332)] # finetune
    q_save = [1.0, 1.0]
    p_save = [1.0, 1.0]
    plt.figure(figsize=(8, 6))
    plt.title('Henon and Heiles  figure 4 var')
    plt.xlabel('y')
    plt.ylabel('vy')
    plt.grid(False)
    # Set color cycle to Tableau 10
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab10.colors)
    for y, vy in y_vy:
        q, p = init_cond_henonheiles(y, vy, H0)

        xplt = []
        yplt = []
        t = 0.0
        while t<tend:
            #print("%e %.4e" % (t, calc_E(q, p)-H0))
            q_save = q.copy()
            p_save = p.copy()
            #kick_inplace(q, p, dt)
            #drift_inplace(q, p, dt)
            triple_jump_inplace(q, p, dt)
            #DKD_LF_inplace(q, p, dt)
            #print("%12.4f %8.3e %8.3e %8.3e %8.3e" % (t, q[0], p[0], q[1], p[1]))
            #print("     %8.3e %8.3e %8.3e" % (q_prev[1], q[1], q_prev[1]*q[1]))
            if q_save[0]*q[0] <= 0 and p[0]>0: # x passed through zero
                xplt.append((q_save[1]+q[1])/2) # we can be more fancy here
                                                # and make linear interpolation,
                                                #let's postpone though
                yplt.append((p_save[1]+p[1])/2)
                print("%12.4f %8.3e %8.3e %.4e" % (t, xplt[-1], yplt[-1],
                                              (calc_E(q, p)-H0)), flush=True)
            t += dt

        plt.scatter(xplt, yplt, s=2, alpha=0.5)
    plt.savefig("my_plot.png")
    plt.show()

if __name__ == "__main__":
    main()
