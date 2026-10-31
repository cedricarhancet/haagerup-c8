"""General exact obstruction-principle machinery (shared by the weight-two and weight-four computations).

W: holomorphic form of weight k for rho = rho_8 with W_z in q^{mu_z} C[[q]], mu_z = Lam_z - 1 (z in P)
   or Lam_z (z not in P).  A 'position' (z, m) denotes the coefficient [q^{Lam_z + m}] W_z.
G: weakly holomorphic form of weight 2 - k for rho^*, G_z = H_{pi z}, H = eta^{-8} Y, Y of weight 6 - k
   for rho.  Then sum_z W_z G_z has weight two and constant term zero (residue theorem), i.e.
       sum_z sum_m W_z[m] G_z[-k_z - m] = 0,   where Lam_{pi z} - 1/3 = -Lam_z + k_z.
Given a finite set A of allowed positions, we find all G (in the span of the given generators) which
pair trivially with every position outside A, and return the induced relations on A.
"""
from fractions import Fraction as Fr
import vvmf_exact as V

d, Lam, names = V.d, V.Lam, V.names

def Y_generators(weight, MJ):
    """generators of weakly holomorphic forms of the given weight (0, 2 or 4) for rho, exact q-series"""
    gens = []
    for j in range(d):
        for m in range(MJ + 1):
            B = V.scal(V.Jpow(m), V.col(j))
            if weight == 0:
                gens.append(('Xi_%s*J^%d' % (names[j], m), B))
                continue
            DB = V.D(B)
            if weight == 2:
                gens.append(('D(Xi_%s*J^%d)' % (names[j], m), DB))
                gens.append(('DJ*Xi_%s*J^%d' % (names[j], m), V.scal(V.DJ, B)))
            elif weight == 4:
                gens.append(('E4*Xi_%s*J^%d' % (names[j], m), V.scal(V.E4, B)))
                gens.append(('serre(D(Xi_%s*J^%d))' % (names[j], m), V.serre(DB, 2)))
                gens.append(('DJ*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.DJ, DB)))
                gens.append(('DJ^2*Xi_%s*J^%d' % (names[j], m), V.scal(V.smul(V.DJ, V.DJ), B)))
            elif weight == 6:
                gens.append(('E6*Xi_%s*J^%d' % (names[j], m), V.scal(V.E6, B)))
                gens.append(('E4*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.E4, DB)))
                S2 = V.serre(DB, 2)
                gens.append(('serre(serre(D(Xi_%s*J^%d)))' % (names[j], m), V.serre(S2, 4)))
                gens.append(('DJ*serre(D(Xi_%s*J^%d))' % (names[j], m), V.scal(V.DJ, S2)))
                gens.append(('DJ*E4*Xi_%s*J^%d' % (names[j], m), V.scal(V.smul(V.DJ, V.E4), B)))
                gens.append(('DJ^2*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.smul(V.DJ, V.DJ), DB)))
                gens.append(('DJ^3*Xi_%s*J^%d' % (names[j], m), V.scal(V.smul(V.smul(V.DJ, V.DJ), V.DJ), B)))
            elif weight == -2:
                gens.append(('E4E6/Delta*Xi_%s*J^%d' % (names[j], m), V.scal(V.E10_over_Delta, B)))
                gens.append(('E4^2/Delta*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.E4sq_over_Delta, DB)))
                gens.append(('E6/Delta*serre(D(Xi_%s*J^%d))' % (names[j], m), V.scal(V.E6_over_Delta, V.serre(DB, 2))))
            elif weight == -4:
                gens.append(('E4^2/Delta*Xi_%s*J^%d' % (names[j], m), V.scal(V.E4sq_over_Delta, B)))
                gens.append(('E6/Delta*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.E6_over_Delta, DB)))
                gens.append(('E4/Delta*serre(D(Xi_%s*J^%d))' % (names[j], m), V.scal(V.E4_over_Delta, V.serre(DB, 2))))
            elif weight == -6:
                gens.append(('E6/Delta*Xi_%s*J^%d' % (names[j], m), V.scal(V.E6_over_Delta, B)))
                gens.append(('E4/Delta*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.E4_over_Delta, DB)))
                gens.append(('1/Delta*serre(serre(D(Xi_%s*J^%d)))' % (names[j], m), V.scal(V.inv_Delta, V.serre(V.serre(DB, 2), 4))))
            elif weight == -8:
                gens.append(('E4/Delta*Xi_%s*J^%d' % (names[j], m), V.scal(V.E4_over_Delta, B)))
                gens.append(('1/Delta*serre(D(Xi_%s*J^%d))' % (names[j], m), V.scal(V.inv_Delta, V.serre(DB, 2))))
            else:
                raise ValueError(weight)
    return gens

def relations(k, P, A, MJ=1, return_G=False):
    """exact relations among the positions A (list of (z, m)) valid for every W in M_k(rho, mu_P)"""
    gens = Y_generators(6 - k, MJ)
    Hs = [V.H_from_Y(Y) for _, Y in gens]
    Aset = set(A)
    kill = []
    for z in range(d):
        kz = V.kz(z)
        mmin = -1 if z in P else 0
        nmin = min([min(H[V.PI[z]].keys()) for H in Hs if H[V.PI[z]]] + [10**6])
        # positions (z, m), m >= mmin, pair with n = -kz - m; we need n >= nmin, i.e. m <= -kz - nmin
        for m in range(mmin, -kz - nmin + 1):
            if (z, m) not in Aset:
                kill.append((z, -kz - m))
    rows = [[H[V.PI[z]].get(n, Fr(0)) for H in Hs] for (z, n) in kill]
    null = V.nullspace(rows, len(Hs))
    rels = []
    for vec in null:
        rels.append([sum(c*H[V.PI[z]].get(-V.kz(z) - m, Fr(0)) for c, H in zip(vec, Hs) if c != 0) for (z, m) in A])
    R, piv = V.rref(rels, len(A)) if rels else ([], [])
    if return_G:
        return R, null, gens
    return R

def exact_forms(k, P, A, MJ=1):
    """exact elements of M_k(rho, mu_P) from the generators; returns their values on the positions A"""
    gens = Y_generators(k, MJ)
    fk = []
    for z in range(d):
        lo = -1 if z in P else 0
        nmin = min([min(Wf[z].keys()) for _, Wf in gens if Wf[z]] + [lo])
        for n in range(nmin, lo):
            fk.append((z, n))
    rows = [[Wf[z].get(n, Fr(0)) for _, Wf in gens] for (z, n) in fk]
    null = V.nullspace(rows, len(gens))
    vals = []
    for vec in null:
        vals.append([sum(c*Wf[z].get(m, Fr(0)) for c, (_, Wf) in zip(vec, gens) if c != 0) for (z, m) in A])
    return vals
