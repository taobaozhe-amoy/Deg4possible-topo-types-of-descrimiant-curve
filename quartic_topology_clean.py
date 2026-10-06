
"""Exact real-projective quartic classifier, usable from SageMath or Python.

Sage usage:
    load('quartic_topology_clean.py')
    answer = classify_delta(A, H, K, L)

Convention:
    q0   = x^T A x + y^T H y
    qinf = 2 x^T L y + 2 y^T K y
    F(p) = (p^T L^T A^(-1) L p)(p^T H p) + (p^T K p)^2.

Rational coefficients only. No plotting, floating-point root matching, CAD
package, or external program is used. The engine uses SymPy's exact polynomial
arithmetic and rational isolating intervals. SymPy can be imported in Sage.

For a reduced quartic with no real singular points, a successful generic
projection certifies one of the six real schemes. Non-real singularities are
allowed and reported separately. Real singularities are NOT forced into the six
smooth schemes. A finite projection-search limit may return NOT_CLASSIFIED.
"""
import sympy as _sp
import random as _random
import builtins as _builtins

_QV = _sp.symbols('u v w')
_X, _Y = _sp.symbols('_X _Y')
_TYPE_TABLE = {
    (0, 1): ('empty', 'empty set', 0),
    (1, 2): ('one_oval', 'one oval', 1),
    (1, 4): ('two_unnested', 'two oval, one is not inside the inteorior of another', 2),
    (2, 3): ('two_nested', 'two ovals ,one in another', 2),
    (1, 6): ('three_ovals', 'three ovals', 3),
    (1, 8): ('four_ovals', 'four ovals', 4),
}


def _rational(a):
    if isinstance(a, (float, _sp.Float)):
        raise ValueError('Please use exact rational coefficients, not floating-point numbers.')
    if hasattr(a, 'parent'):
        par = a.parent()
        if hasattr(par, 'is_exact') and not par.is_exact():
            raise ValueError('Please use QQ matrices, not RR/RDF floating-point matrices.')
    ans = _sp.sympify(str(a))
    if ans.is_Rational is not True:
        raise ValueError('Please use exact rational coefficients, not floating-point numbers.')
    return ans


def _matrix3(M, name):
    if hasattr(M, 'nrows'):  # Sage matrix
        nr, nc = int(M.nrows()), int(M.ncols())
        data = [[M[i, j] for j in range(nc)] for i in range(nr)]
    elif isinstance(M, _sp.MatrixBase):
        nr, nc = M.shape
        data = M.tolist()
    else:
        data = [list(row) for row in M]
        nr = len(data)
        nc = len(data[0]) if nr else 0
    if nr != 3 or nc != 3 or _builtins.any(len(row) != 3 for row in data):
        raise ValueError(name + ' must be a 3×3 matrix.')
    return _sp.Matrix([[_rational(a) for a in row] for row in data])


def _primitive(expr, gens):
    poly = _sp.Poly(expr, *gens, domain=_sp.QQ)
    if poly.is_zero:
        return poly
    _, poly = poly.clear_denoms()
    _, poly = poly.primitive()  # division by a positive content only
    return poly


def _psign(a):
    return 1 if a > 0 else (-1 if a < 0 else 0)


class _AlgebraicSign:
    """Signs at one real root alpha using only exact rational intervals."""
    def __init__(self, squarefree_poly, interval):
        self.P = squarefree_poly
        self.lo, self.hi = map(_sp.Rational, interval)
        if self.lo != self.hi:
            if self.P.eval(self.lo) == 0 or self.P.eval(self.hi) == 0:
                raise ArithmeticError('Root isolation interval endpoints cannot contain another root.')
            if self.P.count_roots(self.lo, self.hi) != 1:
                raise ArithmeticError('Root isolation interval does not exactly isolate one real root.')

    def sign(self, g):
        g = _sp.Poly(g, self.P.gen, domain=_sp.QQ).rem(self.P)
        if g.is_zero:
            return 0
        if self.lo == self.hi:
            return _psign(g.eval(self.lo))
        common = self.P.gcd(g)
        if common.degree() > 0 and common.count_roots(self.lo, self.hi):
            return 0
        # Horner interval evaluation. Refine until zero is excluded.
        coeffs = g.all_coeffs()
        while True:
            low = high = coeffs[0]
            for c in coeffs[1:]:
                candidates = (low*self.lo, low*self.hi,
                              high*self.lo, high*self.hi)
                low, high = _builtins.min(candidates) + c, _builtins.max(candidates) + c
            if low > 0:
                return 1
            if high < 0:
                return -1
            # Exact rational bisection works even when the original
            # isolating interval was produced for a factor of P.
            # P is squarefree and this closed interval has one root.
            for _ in range(4):
                mid = (self.lo + self.hi) / 2
                fm = self.P.eval(mid)
                if fm == 0:
                    self.lo = self.hi = mid
                    return _psign(g.eval(mid))
                if _psign(fm) == _psign(self.P.eval(self.lo)):
                    self.lo = mid
                else:
                    self.hi = mid


class _DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]


def _sphere_regions(counts, events):
    """Glue the complement on S^2, not merely in one affine chart.

    counts[i] = number of real roots in vertical slab i.
    events[i] = (birth, a), where a simple roots lie below the double root.
    sheet 0 represents (X,Y,1), sheet 1 represents its antipode.
    """
    ids, signs = {}, []
    for i, r in enumerate(counts):
        for j in range(r + 1):
            for sheet in (0, 1):
                ids[i, j, sheet] = len(signs)
                signs.append(1 if j % 2 == 0 else -1)
    poles = (len(signs), len(signs) + 1)
    signs += [1, 1]
    dsu = _DSU(len(signs))

    for i, (birth, a) in enumerate(events):
        small = i if birth else i + 1
        big = i + 1 if birth else i
        assert counts[big] == counts[small] + 2
        for j in range(counts[small] + 1):
            targets = [j] if j < a else ([a, a + 2] if j == a else [j + 2])
            for k in targets:
                for sh in (0, 1):
                    dsu.union(ids[small, j, sh], ids[big, k, sh])

    # Infinity: (X,Y,1)/X -> (1,Y/X,0).
    # Passing from +infinity to -infinity reverses root order AND the lift.
    last, r = len(counts) - 1, counts[0]
    assert counts[last] == r
    for j in range(r + 1):
        for sh in (0, 1):
            dsu.union(ids[last, j, sh], ids[0, r-j, 1-sh])

    # Add the two missing points above the projection centre [0:1:0].
    for i, r in enumerate(counts):
        for sh in (0, 1):
            dsu.union(ids[i, r, sh], poles[sh])
            dsu.union(ids[i, 0, sh], poles[1-sh])

    labels = {}
    for j, sg in enumerate(signs):
        root = dsu.find(j)
        if root in labels and labels[root] != sg:
            raise ArithmeticError('Sign region gluing consistency check failed.')
        labels[root] = sg
    return (_builtins.sum(sg == 1 for sg in labels.values()),
            _builtins.sum(sg == -1 for sg in labels.values()))


def _complex_smooth(F):
    """Exact projective smoothness test in all three affine charts."""
    u, v, w = _QV
    derivs = [_sp.diff(F, a) for a in _QV]
    for coord in _QV:
        others = [a for a in _QV if a != coord]
        G = _sp.groebner([g.subs(coord, 1) for g in derivs],
                         *others, order='grevlex', domain=_sp.QQ)
        if not _builtins.any(p.total_degree() == 0 and not p.is_zero for p in G.polys):
            return False
    return True


def _candidate(F, T):
    """Return a certified projection or None if this projection is unsuitable."""
    u, v, w = _QV
    cols = T * _sp.Matrix([u, v, w])
    expr = F.subs(dict(zip(_QV, cols)), simultaneous=True)
    P = _primitive(expr, _QV)
    c = P.coeff_monomial(v**4)
    if c == 0:
        return None
    if c < 0:
        P = -P
        c = -c
    expr = P.as_expr()
    finf = _sp.Poly(expr.subs({u: 1, v: _Y, w: 0}), _Y, domain=_sp.QQ)
    if finf.degree() != 4 or finf.gcd(finf.diff()).degree() > 0:
        return None
    f = _sp.Poly(expr.subs({u: _X, v: _Y, w: 1}), _Y,
                 domain=_sp.QQ.poly_ring(_X))
    D = _primitive(f.discriminant(), (_X,))
    if D.is_zero:
        return None
    # The list interface with strict=True prevents an isolating interval
    # from sharing an endpoint with a different rational root.
    raw_intervals = _sp.intervals([D.as_expr()], eps=_sp.Rational(1, 8), strict=True)
    intervals = [(bounds, multiplicities[0]) for bounds, multiplicities in raw_intervals]
    R = D.sqf_part().set_domain(_sp.QQ)
    subres = f.subresultants(f.diff())
    linear = [g for g in subres if g.degree() == 1]
    if intervals and not linear:
        return None
    if not intervals:
        n = int(_sp.Poly(f.as_expr().subs(_X, 0), _Y).count_roots())
        ni = int(finf.count_roots())
        if n != ni or n not in (0, 2, 4):
            raise ArithmeticError('Critical value fiber count check failed.')
        return {'counts': [n], 'events': [], 'critical_intervals': [],
                'projection': T, 'projected_F': expr, 'projection_discriminant': D.as_expr()}

    lin = linear[-1]
    U = _sp.Poly(lin.nth(1), _X, domain=_sp.QQ)
    N = _sp.Poly(-lin.nth(0), _X, domain=_sp.QQ)
    # beta = N(alpha)/U(alpha) at a simple double-root fibre.
    coeff = [_sp.Poly(f.nth(i), _X, domain=_sp.QQ) for i in range(5)]
    a2, a3 = coeff[2], coeff[3]
    q1n = a3*U + 2*c*N
    q0n = a2*U**2 + 2*a3*N*U + 3*c*N**2
    discq = q1n**2 - 4*c*q0n
    qbeta = a2*U**2 + 3*a3*N*U + 6*c*N**2
    qprime = a3*U + 4*c*N
    fxnum = _builtins.sum((coeff[i].diff()*N**i*U**(3-i) for i in range(4)),
                _sp.Poly(0, _X, domain=_sp.QQ))
    polys = [g.rem(R) for g in (U, discq, qbeta, qprime, fxnum)]
    events = []
    for (interval, multiplicity) in intervals:
        sg = _AlgebraicSign(R, interval)
        su = sg.sign(polys[0])
        if su == 0:
            return None  # two tangencies, flex, or a nongeneric projection
        sfx = sg.sign(polys[4]) * su  # denominator U^3
        if sfx == 0:
            return {'real_singular': True, 'projection': T,
                    'projected_F': expr, 'alpha_interval': (sg.lo, sg.hi),
                    'alpha_polynomial': R.as_expr(),
                    'beta_numerator': N.as_expr(), 'beta_denominator': U.as_expr()}
        if multiplicity != 1:
            return None
        sd = sg.sign(polys[1])  # denominator U^2
        sq = sg.sign(polys[2])  # f_YY/2, denominator U^2
        if sd == 0 or sq == 0:
            return None
        if sd < 0:
            a, other_real = 0, 0
        else:
            other_real = 2
            if sq < 0:
                a = 1
            else:
                sp = sg.sign(polys[3]) * su
                if sp == 0:
                    return None
                a = 0 if sp < 0 else 2
        birth = (-sfx * sq) > 0
        events.append((birth, a, other_real))

    x0 = _sp.floor(intervals[0][0][0]) - 1
    n = int(_sp.Poly(f.as_expr().subs(_X, x0), _Y).count_roots())
    counts = [n]
    for birth, a, other_real in events:
        expected = other_real if birth else other_real + 2
        if n != expected:
            raise ArithmeticError('Critical fiber count check failed.')
        n = n + 2 if birth else n - 2
        if n not in (0, 2, 4):
            raise ArithmeticError('Quartic fiber real root count check failed.')
        counts.append(n)
    if counts[-1] != counts[0] or counts[0] != int(finf.count_roots()):
        raise ArithmeticError('Infinite fiber gluing consistency check failed.')
    return {'counts': counts, 'events': [(b, a) for b, a, r in events],
            'critical_intervals': [z[0] for z in intervals], 'projection': T,
            'projected_F': expr, 'projection_discriminant': D.as_expr()}


def classify_quartic(F, max_tries=60, seed=2026, verbose=True,
                     check_complex=True):
    """Classify a homogeneous QQ quartic in u,v,w (six smooth real schemes).

    F may be a SymPy expression, a Sage polynomial, or an explicit polynomial
    string. Use ** or ^ in a string. Never pass untrusted executable strings.
    Returns a dictionary. 'status' is 'ok', 'real_singular', 'nonreduced',
    'not_classified', or 'invalid_polynomial'. Singular schemes are not labelled
    as smooth schemes. 'ok' is a certified result, not a plotting estimate.
    """
    loc = dict(zip(('u', 'v', 'w'), _QV))
    if not isinstance(F, _sp.Basic):
        F = _sp.sympify(str(F).replace('^', '**'), locals=loc)
    F = _sp.expand(F)
    if F.has(_sp.Float):
        raise ValueError('Please use exact rational coefficients, not floating-point numbers.')
    P = _sp.Poly(F, *_QV, domain=_sp.QQ)
    out = {'F': F}
    if P.is_zero or set(_builtins.sum(m) for m, c in P.terms()) != {4}:
        out.update(status='invalid_polynomial', type=' not a non-zero homogeneous quartic')
        return _report(out, verbose)
    common = P
    for t in _QV:
        common = common.gcd(P.diff(t))
    if common.total_degree() > 0:
        out.update(status='nonreduced', type='nonreduced quartic curve: cannot be applied to the six types table')
        return _report(out, verbose)
    out['complex_smooth'] = _complex_smooth(F) if check_complex else None
    rng = _random.Random(int(seed))
    for attempt in range(int(max_tries)):
        if attempt == 0:
            T = _sp.eye(3)
        else:
            bound = 2 + attempt // 20
            while True:
                T = _sp.Matrix(3, 3, [rng.randint(-bound, bound) for _ in range(9)])
                if T.det() != 0:
                    break
        data = _candidate(F, T)
        if data is None:
            continue
        if data.get('real_singular'):
            out.update(data)
            out.update(status='real_singular', type='contains real singularities: does not belong to the six types of smooth real images')
            return _report(out, verbose)
        npos, nneg = _sphere_regions(data['counts'], data['events'])
        key = tuple(sorted((npos, nneg)))
        if key not in _TYPE_TABLE:
            raise ArithmeticError('The number of sphere regions does not match the classification of smooth real quartic curves. Please retain the input for inspection.')
        code, label, ovals = _TYPE_TABLE[key]
        out.update(data)
        out.update(status='ok', type=label, type_code=code, ovals=ovals,
                   sphere_region_counts=(npos, nneg),
                   projection_attempts=attempt + 1,
                   real_smooth=True, exact=True)
        return _report(out, verbose)
    out.update(status='not_classified',
               type='not classified: please increase max_tries; there may also be higher-order real singularities')
    return _report(out, verbose)


def _report(out, verbose):
    if not verbose:
        return out
    print('F =', out['F'])
    if out.get('complex_smooth') is not None:
        print('Delta in C is smooth:', 'Yes' if out['complex_smooth'] else 'No')
    if out['status'] == 'ok':
        print('Sphere region counts (normalized equation):', out['sphere_region_counts'])
        print('Number of real critical values:', len(out['events']))
        print('Classification type:', out['type'])
        print('Number of connected components of Delta(R):', out['ovals'])
        print('Status: Exact Real Root Isolation and Projective Gluing Judgment')
    else:
        print('Classification Result:', out['type'])
    return out


def classify_delta(A, H, K, L, max_tries=60, seed=2026, verbose=True,
                   check_model=True):
    """Sage entry point: input A,H,K symmetric, A invertible, arbitrary L.

    The entries must be exact rationals. By default check that the associated
    intersection of two quadrics is a smooth projective threefold.
    """
    A, H, K, L = [_matrix3(M, name) for M, name in
                  zip((A, H, K, L), ('A', 'H', 'K', 'L'))]
    for name, M in (('A', A), ('H', H), ('K', K)):
        if M != M.T:
            raise ValueError(name + ' must be symmetric.')
    if A.det() == 0:
        raise ValueError('A must be invertible.')
    p = _sp.Matrix(_QV)
    S = (p.T * L.T * A.inv() * L * p)[0]
    h = (p.T * H * p)[0]
    k = (p.T * K * p)[0]
    F = _sp.expand(S*h + k*k)
    model_ok, pencil = None, None
    if check_model:
        t = _sp.Symbol('_t')
        M0 = _sp.diag(A, H)
        Mi = _sp.zeros(6)
        Mi[:3, 3:] = L
        Mi[3:, :3] = L.T
        Mi[3:, 3:] = 2*K
        pencil = _sp.Poly((t*Mi - M0).det(method='domain-ge'), t, domain=_sp.QQ)
        model_ok = (not pencil.is_zero and pencil.degree() in (5, 6)
                    and pencil.gcd(pencil.diff()).degree() == 0)
        if verbose:
            print('Z smoothness check:', 'Passed' if model_ok else 'Failed')
        if not model_ok:
            out = {'F': F, 'status': 'invalid_model', 'model_smooth': False,
                   'type': 'the corresponding Z is not a smooth threefold; not classified as a conic bundle'}
            return _report(out, verbose)
    out = classify_quartic(F, max_tries=max_tries, seed=seed, verbose=verbose)
    out.update(S=_sp.expand(S), h=_sp.expand(h), k=_sp.expand(k),
               model_smooth=model_ok,
               pencil_polynomial=pencil.as_expr() if pencil is not None else None)
    return out
