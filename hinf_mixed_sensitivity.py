"""
Unified Toker-Ozbay / Foias-Ozbay-Tannenbaum H-infinity mixed-sensitivity solver.

Given:
    W1(s) = nw1(s)/dw1(s)      performance weight on S
    W2(s) = nw2(s)/dw2(s)      performance weight on T
    Md(s) = nmd(s)/dmd(s)      finite Blaschke product for the plant's unstable
                                poles (alpha = roots(nmd), the plant's RHP poles)
    Mn(s)                      inner (all-pass) part of the plant not captured by
                                Md -- typically carries the delays. Function handle.
    No(s)                      outer (stable, minimum-phase) part of the plant.
                                Function handle (may also carry delays).

    Plant: P(s) = Mn(s)*No(s)/Md(s)

Returns:
    gamma_opt : optimal H-infinity mixed-sensitivity level
    omg, C, S, T, per : optimal controller's frequency response and closed-loop
                         performance, evaluated on a log-spaced grid.
"""
import numpy as np


# ---------------- generic polynomial helpers (MATLAB-style, descending powers) ----------------
def star(p):
    p = np.asarray(p, dtype=complex)
    n = len(p)
    k = np.arange(1, n+1)
    return p*(-1.0)**(n+k)

def polyadd(p, q):
    n = max(len(p), len(q))
    pp = np.concatenate([np.zeros(n-len(p)), p])
    qq = np.concatenate([np.zeros(n-len(q)), q])
    s = pp+qq
    EPS = 1e-20
    while len(s) > 1 and abs(s[0]) < EPS:
        s = s[1:]
    return s

def conv(a, b):
    return np.convolve(a, b)

def pspec(p):
    """Return the (real-coefficient) polynomial built from the strictly-left-half-plane
    roots of p, plus one root from each imaginary-axis conjugate pair."""
    EPS = 1e-6
    z = np.roots(p)
    zn = z[np.real(z) < -EPS]
    zi = list(z[np.abs(np.real(z)) <= EPS])
    zx = []
    while len(zi) > 0:
        zx.append(zi[0])
        if len(zi) > 1:
            dists = [abs(zi[j]-zi[0]) for j in range(1, len(zi))]
            n = int(np.argmin(dists)) + 1
            zi = zi[1:n] + zi[n+1:]
        else:
            zi = []
    return np.poly(list(zn)+zx)

def spec(nums, dens):
    """Spectral factorization: given a para-Hermitian A(s)/B(s) >= 0 on jw, return
    (dG, nG) such that dG(s)*dG(-s)/(nG(s)*nG(-s)) == A(s)/B(s)/... (see usage in
    build_gamma_matrix). Order matches MATLAB's `[dG,nG] = spec(A,B)` convention:
    dG is the SCALED factor of the first argument (A), nG is the plain factor of the
    second argument (B)."""
    numK = nums[0]/dens[0]*(-1)**((len(nums)-len(dens))/2)
    if numK < 0:
        raise ValueError('Cannot do spectral factorization (numK<0) at this gamma.')
    numK = np.sqrt(numK)
    dG = numK*pspec(nums)
    nG = pspec(dens)
    return np.real(dG), np.real(nG)


class HinfMixedSens:
    def __init__(self, nw1, dw1, nw2, dw2, nmd, dmd, mn, no, EPS=1e-6):
        if not callable(mn):
            raise TypeError(
                "mn must be a callable, e.g. mn = lambda s: np.exp(-0.2*s)*(...)/(...). "
                f"Got {type(mn)} instead -- if you evaluated it at specific s-values and "
                "passed the resulting array, that's the bug: pass the function itself.")
        if not callable(no):
            raise TypeError(
                f"no must be a callable (function of s). Got {type(no)} instead.")
        self.nw1, self.dw1 = np.asarray(nw1, float), np.asarray(dw1, float)
        self.nw2, self.dw2 = np.asarray(nw2, float), np.asarray(dw2, float)
        self.nmd, self.dmd = np.asarray(nmd, float), np.asarray(dmd, float)
        self.mn = mn
        self.no = no
        self.md = lambda s: np.polyval(self.nmd, s)/np.polyval(self.dmd, s)
        self.alpha = np.roots(self.nmd)          # plant's unstable (RHP) poles
        self.EPS = EPS

    # ---- build the Pick matrix at a given gamma; returns (M, nE,dE,nF,dF) or None if non-generic
    def _pick_matrix(self, gamma):
        nw1, dw1, nw2, dw2 = self.nw1, self.dw1, self.nw2, self.dw2
        a1 = polyadd(conv(nw1, star(nw1)), -gamma**2*conv(dw1, star(dw1)))
        b1 = gamma**2*conv(dw1, star(dw1))
        a2 = polyadd(conv(nw2, star(nw2)), -gamma**2*conv(dw2, star(dw2)))
        b2 = gamma**2*conv(dw2, star(dw2))
        A1, B1 = conv(a1, a2), conv(b1, b2)
        A = polyadd(B1, -A1)
        B = B1

        try:
            dG, nG = spec(A, B)
        except ValueError:
            return None

        eta = np.roots(dw1)
        n1 = len(eta)
        nF = conv(nG, np.poly(-eta))
        dF = conv(dG, np.poly(eta))

        beta = np.roots(a1)
        beta = beta[(np.real(beta) > self.EPS) |
                    ((np.abs(np.real(beta)) <= self.EPS) & (np.imag(beta) > 0))]
        alpha = self.alpha
        l = len(alpha)
        if len(beta) != n1:
            return None

        N = n1+l
        rows = []
        for pt in beta:
            powers = pt**np.arange(N-1, -1, -1)
            val = self.mn(pt)*np.polyval(nF, pt)/np.polyval(dF, pt)
            rows.append(np.concatenate([powers, val*powers]))
        for pt in alpha:
            powers = pt**np.arange(N-1, -1, -1)
            val = self.mn(pt)*np.polyval(nF, pt)/np.polyval(dF, pt)
            rows.append(np.concatenate([powers, val*powers]))
        for pt in beta:
            powers = (-pt)**np.arange(N-1, -1, -1)
            val = self.mn(pt)*np.polyval(nF, pt)/np.polyval(dF, pt)
            rows.append(np.concatenate([val*powers, powers]))
        for pt in alpha:
            powers = (-pt)**np.arange(N-1, -1, -1)
            val = self.mn(pt)*np.polyval(nF, pt)/np.polyval(dF, pt)
            rows.append(np.concatenate([val*powers, powers]))
        M = np.array(rows, dtype=complex)
        return dict(M=M, nF=nF, dF=dF, a1=a1, b1=b1, n1=n1, l=l, N=N)

    def smin(self, gamma):
        sol = self._pick_matrix(gamma)
        if sol is None:
            return None
        s = np.linalg.svd(sol['M'], compute_uv=False)
        return s[-1]

    # ---- automated search for gamma_opt: the LARGEST gamma at which the Pick
    # matrix is (numerically) singular.
    def find_gamma_opt(self, gmin=1e-3, gmax=20.0, n_coarse=4000, refine_tol=1e-10):
        gam = np.linspace(gmin, gmax, n_coarse)
        vals = np.full(n_coarse, np.nan)
        for i, g in enumerate(gam):
            v = self.smin(g)
            if v is not None:
                vals[i] = v

        # find all local minima among the generic (non-nan) samples
        candidates = []
        for i in range(1, n_coarse-1):
            if np.isnan(vals[i]):
                continue
            left = vals[i-1] if not np.isnan(vals[i-1]) else np.inf
            right = vals[i+1] if not np.isnan(vals[i+1]) else np.inf
            if vals[i] <= left and vals[i] <= right and vals[i] < 1e-2:
                candidates.append(gam[i])

        if not candidates:
            raise RuntimeError('No near-singular gamma found in the scanned range; '
                                'widen [gmin,gmax] or check the plant/weight data.')

        # refine each candidate with a local golden-section-style bisection
        refined = []
        for g0 in candidates:
            lo, hi = max(gmin, g0*0.98), min(gmax, g0*1.02)
            for _ in range(60):
                mid1 = lo + (hi-lo)/3
                mid2 = hi - (hi-lo)/3
                v1, v2 = self.smin(mid1), self.smin(mid2)
                if v1 is None: v1 = np.inf
                if v2 is None: v2 = np.inf
                if v1 < v2:
                    hi = mid2
                else:
                    lo = mid1
                if hi-lo < refine_tol:
                    break
            gstar = (lo+hi)/2
            refined.append(gstar)

        gamma_opt = max(refined)   # per Toker-Ozbay: take the LARGEST singular gamma
        return gamma_opt, sorted(refined)

    # ---- build optimal controller + frequency response at a given (already found) gamma
    def controller_response(self, gamma, lw1=-4, lw2=4, Lpoints=2500):
        nw1, dw1, nw2, dw2 = self.nw1, self.dw1, self.nw2, self.dw2

        # gamma sits exactly on the boundary where the Pick matrix becomes singular;
        # at that exact floating-point value root-finding can occasionally flip to a
        # non-generic root count (especially when there are few/no unstable poles).
        # If that happens, nudge gamma by a tiny amount until a generic point is found.
        sol = self._pick_matrix(gamma)
        if sol is None:
            nudged = False
            for delta in [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]:
                for sgn in (1, -1):
                    g_try = gamma + sgn*delta*max(1, gamma)
                    sol_try = self._pick_matrix(g_try)
                    if sol_try is not None:
                        print(f'hinf_mixed_sensitivity: gamma_opt={gamma:.10g} was exactly on a '
                              'non-generic boundary (common when the plant has few/no unstable '
                              f'poles); nudged to {g_try:.10g} to build the controller.')
                        gamma = g_try
                        sol = sol_try
                        nudged = True
                        break
                if nudged:
                    break
            if sol is None:
                raise RuntimeError(
                    f'Pick matrix is non-generic at and near gamma={gamma} and could not be '
                    'nudged to a generic point. Try a finer n_coarse, or a tighter [gmin,gmax] '
                    'bracket around this value.')

        a1 = polyadd(conv(nw1, star(nw1)), -gamma**2*conv(dw1, star(dw1)))
        b1 = gamma**2*conv(dw1, star(dw1))
        nE, dE = a1, b1

        M, nF, dF, N = sol['M'], sol['nF'], sol['dF'], sol['N']

        Muu, Suu, Vuuh = np.linalg.svd(M)
        Vuu = Vuuh.conj().T
        LL = Vuu[:, -1]
        smin = Suu[-1]
        nL = np.real(LL[N:2*N])
        dL = np.real(LL[0:N])

        lws = np.linspace(lw1, lw2, Lpoints+1)
        omg = 10.0**lws
        s = 1j*omg
        c1 = np.polyval(nE, s)/np.polyval(dE, s)
        c2 = np.polyval(nF, s)/np.polyval(dF, s)
        c3 = np.polyval(nL, s)/np.polyval(dL, s)
        C = c1*c2*c3/(1+self.mn(s)*c2*c3)*self.md(s)/self.no(s)

        P = self.mn(s)*self.no(s)/self.md(s)
        S = 1/(1+P*C)
        T = 1-S
        W1v = np.polyval(nw1, s)/np.polyval(dw1, s)
        W2v = np.polyval(nw2, s)/np.polyval(dw2, s)
        per = np.sqrt(np.abs(W1v*S)**2 + np.abs(W2v*T)**2)

        return dict(omg=omg, C=C, S=S, T=T, per=per, smin=smin, gamma=gamma,
                    nE=nE, dE=dE, nF=nF, dF=dF, nL=nL, dL=dL)


def solve_hinf_mixed_sensitivity(nw1, dw1, nw2, dw2, nmd, dmd, mn, no,
                                  gmin=1e-3, gmax=20.0, n_coarse=4000):
    """Top-level convenience function: returns (gamma_opt, omg, C, S, T, per, candidates).

    ALWAYS sanity-check the returned `per` array before trusting gamma_opt: for a
    genuinely optimal solution it must be FLAT and equal to gamma_opt at every
    frequency (the "all-pass completion" signature of true H-infinity optimality).
    This function checks that automatically and raises a warning if it isn't.
    """
    prob = HinfMixedSens(nw1, dw1, nw2, dw2, nmd, dmd, mn, no)
    gamma_opt, candidates = prob.find_gamma_opt(gmin=gmin, gmax=gmax, n_coarse=n_coarse)
    resp = prob.controller_response(gamma_opt)
    gamma_opt = resp['gamma']   # may have been nudged off a non-generic boundary point
    per = resp['per']
    if (per.max() - per.min()) > 1e-3*gamma_opt:
        import warnings
        warnings.warn(
            f"per is not flat (min={per.min():.6g}, max={per.max():.6g}) -- this "
            "gamma_opt may not be the true optimum. Check gmin/gmax bracket and EPS.")
    return gamma_opt, resp['omg'], resp['C'], resp['S'], resp['T'], per, candidates


if __name__ == '__main__':
    # ----------------------------------------------------------------------
    # Validated example (Ozbay, Gumussoy, Kashima, Yamamoto, "Frequency Domain
    # Techniques for H-infinity Control of Distributed Parameter Systems",
    # SIAM 2018, Section 6.3.1, Eq. 6.62): unstable retarded time-delay plant
    #     P(s) = (s-1)(s+4)e^{-hs} / [(s^2+8s+17)(s+1-3e^{-2hs})],  h=ln(2)
    # Weights: W1(s)=1/s, W2(s)=2s.  Book reports gamma_opt = 17.846.
    # ----------------------------------------------------------------------
    h = np.log(2)
    nw1 = [1.0];        dw1 = [1.0, 0.0]      # W1 = 1/s
    nw2 = [2.0, 0.0];   dw2 = [1.0]           # W2 = 2*s
    nmd = [1.0, -0.5];  dmd = [1.0, 0.5]      # Md(s) = (s-0.5)/(s+0.5)

    def mn(s):
        return (s-1)/(s+1)*np.exp(-h*s)

    def no(s):
        return (2*(s+4)*(s+1)/((s**2+8*s+17)*(s+0.5))
                * (s-0.5)/(s+1-3*np.exp(-2*h*s)))

    gamma_opt, omg, C, S, T, per, candidates = solve_hinf_mixed_sensitivity(
        nw1, dw1, nw2, dw2, nmd, dmd, mn, no, gmin=1.0, gmax=30.0, n_coarse=3000)

    print('candidate near-singular gammas:', candidates)
    print(f'gamma_opt = {gamma_opt:.6f}   (book reports 17.846)')
    print(f'per: min={per.min():.6f}  max={per.max():.6f}  (should be flat)')

