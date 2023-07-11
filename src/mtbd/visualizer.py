import numpy as np
from bdpn.mtbd import find_index_within_bounds
from matplotlib import pyplot as plt
from scipy.integrate import odeint



def plot_P_simple(k, get_U, ti, t0, MU, LA, PSI):
    t = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    plt.plot(t, sol[:, 0], 'b', label='logP_E(t)')
    plt.plot(t, sol[:, 1], 'g', label='logP_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()



def plot_U(get_U, T):
    t = np.linspace(0, T, 1001)
    esol_E = [get_U(_)[0] for _ in t]
    esol_I = [get_U(_)[1] for _ in t]
    plt.plot(t, esol_E, 'b', label='U_E(t)')
    plt.plot(t, esol_I, 'g', label='U_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def plot_P(k, get_U, ti, t0, MU, LA, PSI):
    tt = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, tt)

    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI

    def get_P_Euler():

        def pdf_Pany_l(P, t):
            U = get_U(t)
            return (SIGMA - LA.dot(U)) * P - (MU + U * LA).dot(P)

        dt = 1e-3
        tj = ti
        yj = np.array(y0)
        cs = [1]
        ts = [ti]
        ys = [y0]

        while tj > t0:
            yj_next = yj - dt * pdf_Pany_l(yj, tj)
            if np.any(yj_next < 0) or np.all(yj_next == 0) or np.any(yj_next > np.prod(cs)):
                c = 1 / min(yj[yj > 0])
                cs.append(c)
                yj *= c
                continue
            else:
                tj -= dt
                yj = yj_next
                ts.append(tj)
                ys.append(yj / np.prod(cs))
        print(cs)
        return np.array(ts), np.array(ys)

    tjs, sol3 = get_P_Euler()

    i = find_index_within_bounds(sol, 0, len(tt), 1)
    vs = sol[i, :]
    c = 1 / min(sol[i, sol[i, :] > 0])
    print(c)
    y0 = vs * c
    ti = tt[i]
    nsteps = 100
    tt2 = np.linspace(ti, t0, nsteps)
    sol2 = odeint(pdf_Pany_l, y0, tt2) / c

    plt.plot(tt, sol[:, 0], 'b', label='P_E(t)')
    plt.plot(tt, sol[:, 1], 'g', label='P_I(t)')
    plt.plot(tt2, sol2[:, 0], 'bv', label='P_E(t)')
    plt.plot(tt2, sol2[:, 1], 'gv', label='P_I(t)')
    plt.plot(tjs, sol3[:, 0], 'b*', label='P_E(t)')
    plt.plot(tjs, sol3[:, 1], 'g*', label='P_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()