import numpy as np

import numba



@numba.jit
def calculate_a_soave(T: float, va0, vc1, vtc):
    vtr = T / vtc[:]
    va = va0[:] * (1 + vc1[:] * (1 - np.sqrt(vtr[:])))**2
    return va

@numba.jit
def calculate_a_mix(vx, T: float, rho: float, va, kij, ncomp):
    V = 1 / rho
    am = 0.
    for i in range(ncomp):
        for j in range(ncomp):
            aij = np.sqrt(va[i] * va[j]) * (1 - kij[i, j])   # Parâmetro cruzado
            am += vx[i] * vx[j] * aij
    return am

@numba.jit
def calculate_b_mix(vx, vb):
    return np.sum(vx * vb)
@numba.jit
def calculate_dadn(vx, T: float, rho: float, va, am, kij, ncomp):
    dadn = np.zeros(ncomp)
    sum1 = np.zeros(ncomp)

    for i in range(ncomp):
        for j in range(ncomp):
            aij = np.sqrt(va[i] * va[j]) * (1 - kij[i, j])
            sum1[i] += vx[j] * aij
        dadn[i] = 2 * sum1[i] - am
    return dadn
@numba.jit
def calculate_pressure(vx, T: float, rho: float, am, bm, _R, sig, eps):
    P = (_R * T) / (1 / rho - bm) - am / ((1 / rho + sig * bm) * (1 / rho + eps * bm))
    return P
@numba.jit
def calculate_phi(vx, T: float, rho: float, am, bm, da, db, P, _R, sig, eps, ncomp):
    V = 1 / rho
    lnPhi = np.zeros(ncomp)
    I = (1. / (eps - sig)) * np.log((V + eps * bm) / (V + sig * bm))
    q = am / (bm * _R * T)
    dq = q * (1 + da[:]/am - db[:]/bm)

    for i in range(ncomp):
        lnPhi[i] = (
            (db[i] / bm) * ((P * V) / (_R * T) - 1.) -
            np.log(P * (V - bm) / (_R * T)) -
            dq[i] * I
        )

    return np.exp(lnPhi)
@numba.jit
def calculate_G_res_partial(vphi, _R, T):
    vlnphi = np.log(vphi)
    vGres = (_R * T) * vlnphi
    return vGres
@numba.jit
def calculate_G_res_mix(vGres, vx):
    return np.sum(vGres * vx)
