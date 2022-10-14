import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel, BirthDeathModel
from mtbd.mtbd_estimator import find_index_within_bounds, state_frequencies, compute_U


def plot_P_simple(k, get_U, ti, t0, MU, LA, PSI):
    t = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    plt.plot(t, sol[:, 0], 'b', label='logP(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def plot_U(get_U, T):
    t = np.linspace(0, T, 1001)
    esol_E = [get_U(_)[0] for _ in t]
    plt.plot(t, esol_E, 'b', label='U(t)')
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
    y0[0] = 1
    sol = odeint(pdf_Pany_l, y0, tt)


    i = find_index_within_bounds(sol, 0, len(tt), 1)
    vs = sol[i, :]
    c = 1 / min(sol[i, sol[i, :] > 0])
    print(c)
    y0 = vs * c
    ti = tt[i]
    nsteps = 100
    tt2 = np.linspace(ti, t0, nsteps)
    sol2 = odeint(pdf_Pany_l, y0, tt2) / c


    plt.plot(tt, sol[:, 0], 'b', label='logP_E(t)')
    plt.plot(tt2, sol2[:, 0], 'bv', label='logP_E(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


if __name__ == '__main__':

    p = 0.6225767763228239
    real_la = 0.907081384137969
    real_psi = 0.2692907505391973
    T = 20
    model = BirthDeathModel(la=real_la, psi=real_psi, p=p)

    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    PI = state_frequencies(MU, LA, PSI)
    PSI_RHO = PSI * RHO
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)

    plot_U(get_U, T)

    plot_P(1, get_U, ti=10, t0=5, MU=MU, LA=LA, PSI=PSI)