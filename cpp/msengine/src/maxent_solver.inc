// RESEARCH (MS_SOLVER, MS_SOLVER_CMP): other solvers for the fixed alpha, fixed envelope problem of pass 1
//
//   Q(f) = chi2 / 2 + alpha sum_j (f_j ln(f_j / m) - f_j + m),   f = exp(u), u in [lo, hi]
//
// the problem lbfgs_b solves (Objective: chi2 / 2 - alpha x entropy). Included into maxent.cpp after lbfgs_b.
//   mm       majorise-minimise: with H = A'WA, c = A'Wy, d = (H f) / f, the separable majoriser of the data term plus
//            the unchanged entropy has the minimiser ln x = ln m + c / alpha - omega(ln(d m / alpha) + c / alpha)
//            (omega: Wright omega); one forward and one adjoint product per step, Q never rises
//   mmacc    mm with momentum (FISTA in u) and a gradient restart
//   squarem  mm accelerated by SQUAREM (Varadhan and Roland 2008, scheme S3) with a fall back to plain steps
//   lbfgsmm  lbfgs_b with the majoriser's curvature as its first inverse Hessian (1 / (f (H f + alpha)) in u)
//   tn       truncated Newton: (H + alpha / f) df = -grad by conjugate gradients preconditioned with the majoriser's
//            diagonal (d + alpha / f), the step in u = df / f, backtracking on Q
//   dual     Newton-CG on the dual: lambda in data space, f(lambda) = m exp(-A'lambda / alpha) (clipped to the box)
//
// The duality gap (Fenchel) of a primal f with lambda = W (A f - y) is computed for every result: it bounds Q(f) - Q*.

inline double flog(double x) {   // fdlibm e_log.c (error below 1 ulp) for normal positive x; std::log otherwise
    if (!(x >= 2.2250738585072014e-308) || !(x < 1.7976931348623157e308)) return std::log(x);
    const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10,
                 Lg1 = 6.666666666666735130e-01, Lg2 = 3.999999999940941908e-01, Lg3 = 2.857142874366239149e-01,
                 Lg4 = 2.222219843214978396e-01, Lg5 = 1.818357216161805012e-01, Lg6 = 1.531383769920937332e-01,
                 Lg7 = 1.479819860511658591e-01;
    uint64_t b;
    std::memcpy(&b, &x, 8);
    int32_t hx = (int32_t)(b >> 32);
    int k = (hx >> 20) - 1023;
    hx &= 0x000fffff;
    const int32_t i0 = (hx + 0x95f64) & 0x100000;
    b = (b & 0xffffffffull) | ((uint64_t)(uint32_t)(hx | (i0 ^ 0x3ff00000)) << 32);
    double xx;
    std::memcpy(&xx, &b, 8);
    k += (i0 >> 20);
    const double f = xx - 1.0;
    const double s = f / (2.0 + f), dk = (double)k, z = s * s, w = z * z;
    int32_t i = hx - 0x6147a;
    const int32_t j = 0x6b851 - hx;
    const double t1 = w * (Lg2 + w * (Lg4 + w * Lg6)), t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
    i |= j;
    const double R = t2 + t1;
    if (i > 0) {
        const double hfsq = 0.5 * f * f;
        return dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);
    }
    return dk * ln2_hi - ((s * (f - R) - dk * ln2_lo) - f);
}

// Wright omega: the w > 0 with w + ln w = t; w0 > 0: a start (the previous step's value). From a start within
// |r| < 0.05 one step of Fritsch, Shafer and Crowley (error of order r^4) is enough; else up to four
inline double wright_omega(double t, double w0) {
    if (t < -36.0) return fexp(t);          // w = e^t (1 - e^t + ...): below 1e-15 relative
    double w = w0, r = 1.0;
    if (w > 0) r = t - w - flog(w);
    if (!(w > 0) || std::fabs(r) > 0.5) {
        if (t > 1.0) { const double lt = flog(t); w = t - lt + lt / t; }
        else { const double zz = fexp(t), L = std::log1p(zz); w = L * (1.0 - std::log1p(L) / (2.0 + L)); }
        r = t - w - flog(w);
    }
    for (int it = 0; it < 4; it++) {
        const double q = 2.0 * (1.0 + w) * (1.0 + w + (2.0 / 3.0) * r);
        w = w * (1.0 + r / (1.0 + w) * (q - r) / (q - 2.0 * r));
        if (std::fabs(r) < 0.05) break;
        r = t - w - flog(w);
        if (std::fabs(r) < 1e-15 * (1.0 + std::fabs(t))) break;
    }
    return w;
}

struct ProdCount {   // forward and adjoint products since the start (from the profile counters)
    long f0, a0;
    ProdCount() : f0(prof_n[T_FORWARD]), a0(prof_n[T_ADJOINT]) {}
    long prod() const { return (prof_n[T_FORWARD] - f0) + (prof_n[T_ADJOINT] - a0); }
};

struct SolveStats {
    double Q0 = 0, Q = 0, gap = 0;
    long prod = 0, prod_target = -1;   // products in all, and until Q first reached the target
    int iters = 0;
    double sec = 0;
    double chi2 = 0;
    double Qat[6] = {0, 0, 0, 0, 0, 0};   // Q after 20, 40, 60, 80, 100, 150 products (mmf)
};
static const long kCheck[6] = {20, 40, 60, 80, 100, 150};

struct Inner {
    Model* model; const std::vector<double>* y; const std::vector<double>* w;
    double alpha, dflt, ld, lo, hi;
    size_t na, ng;
    std::vector<double> c;   // A'(w y), for the envelope of the problem
    bool have_c = false;
    Inner(const Objective& o, double lo_, double hi_)
        : model(o.model), y(o.y), w(o.w), alpha(o.alpha), dflt(o.dflt), ld(o.log_dflt), lo(lo_), hi(hi_),
          na(o.model->na), ng(o.y->size()) {}
    void make_c() {
        if (have_c) return;
        std::vector<double> wy(ng);
        const double* yy = y->data(); const double* ww = w->data();
        pfor(ng, [&](size_t a, size_t b) { for (size_t i = a; i < b; i++) wy[i] = ww[i] * yy[i]; });
        c.assign(na, 0.0);
        model->adjoint_f(wy.data(), c.data());
        have_c = true;
    }
    // Q at u (as Objective computes it): f = exp(u) and fit = A f are written
    double value(const std::vector<double>& u, std::vector<double>& f, std::vector<double>& fit) const {
        Timed tm(T_MMVAL);
        f.resize(na); fit.resize(ng);
        double ent;
        const double* uu = u.data(); double* ff = f.data();
        const double dl = dflt, l = ld;
        psums(na, 1, &ent, [&](size_t i0, size_t i1, double* acc) {
            fexp_n(uu, ff, i0, i1);
            double s = 0;
            for (size_t j = i0; j < i1; j++) { const double fj = ff[j]; s += fj - dl - fj * (uu[j] - l); }
            acc[0] = s;
        });
        model->forward(ff, fit.data());
        return chi2half(fit) - alpha * ent;
    }
    double chi2half(const std::vector<double>& fit) const {
        const double* yy = y->data(); const double* ww = w->data(); const double* ft = fit.data();
        return 0.5 * psum(ng, [&](size_t i) { const double d = ft[i] - yy[i]; return ww[i] * d * d; });
    }
    double entropy_term(const std::vector<double>& u, const std::vector<double>& f) const {   // sum f ln(f/m) - f + m
        return psum(na, [&](size_t j) { return f[j] * (u[j] - ld) - f[j] + dflt; });
    }
    // the Fenchel gap of (u, fit = A exp(u)) with lambda = W (fit - y): Q(f) + F(lambda) >= Q(f) - Q*. One adjoint.
    double gap(const std::vector<double>& u, const std::vector<double>& f, const std::vector<double>& fit) const {
        std::vector<double> lam(ng), s(na);
        const double* yy = y->data(); const double* ww = w->data();
        pfor(ng, [&](size_t a, size_t b) { for (size_t i = a; i < b; i++) lam[i] = ww[i] * (fit[i] - yy[i]); });
        model->adjoint_f(lam.data(), s.data());
        const double c2 = chi2half(fit);
        const double Q = c2 + alpha * entropy_term(u, f);
        const double ylam = psum(ng, [&](size_t i) { return yy[i] * lam[i]; });
        const double al = alpha, l = ld, dl = dflt, lo_ = lo, hi_ = hi;
        const double conj = psum(na, [&](size_t j) {
            const double lf = std::min(hi_, std::max(lo_, l - s[j] / al)), fh = fexp(lf);
            return al * (fh * (lf - l) - fh + dl) + s[j] * fh;
        });
        return Q + (c2 + ylam - conj);
    }
};

// ---- majorise-minimise steps
// x = MM(u): h = A'(w A f) at u (given), c = A'(w y); om: Wright omega of the last step (start values)
static void mm_map(const Inner& P, const std::vector<double>& u, const std::vector<double>& h, std::vector<double>& x,
                   std::vector<double>& om) {
    Timed tm(T_MMMAP);
    const double la = std::log(P.alpha), al = P.alpha, ld = P.ld, lo = P.lo, hi = P.hi;
    const double* cc = P.c.data();
    x.resize(P.na);
    if (om.size() != P.na) om.assign(P.na, -1.0);
    pfor(P.na, [&](size_t a, size_t b) {
        for (size_t j = a; j < b; j++) {
            const double ca = cc[j] / al;
            double v;
            if (h[j] > 0) {
                const double t = flog(h[j]) - u[j] + ld - la + ca;
                const double o = wright_omega(t, om[j]);
                om[j] = o;
                v = ld + ca - o;
            } else {
                v = ld + ca;
            }
            x[j] = std::min(hi, std::max(lo, v));
        }
    });
}

// h = A'(w fit)
static void hf_of(const Inner& P, const std::vector<double>& fit, std::vector<double>& h, std::vector<double>& wf) {
    wf.resize(P.ng); h.resize(P.na);
    const double* ww = P.w->data();
    pfor(P.ng, [&](size_t a, size_t b) { for (size_t i = a; i < b; i++) wf[i] = ww[i] * fit[i]; });
    P.model->adjoint_f(wf.data(), h.data());
}

// mode 0: plain MM; 1: momentum with gradient restart
static void mm_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st, int mode, double rtol = 0) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na;
    std::vector<double> f, fit, h, wf, x, om, xprev = u, best;
    double Q = P.value(u, f, fit);
    st.Q0 = Q;
    double bestQ = Q;
    best = u;
    double tk = 1.0;
    int it = 0;
    const double al = P.alpha, ld = P.ld;
    for (; it < maxit && pc.prod() < budget; it++) {
        hf_of(P, fit, h, wf);
        mm_map(P, u, h, x, om);
        if (mode == 1) {
            // restart when the momentum points uphill: grad_u(y) . (x - x_prev) > 0, grad_u = f (h - c + alpha (u - ld))
            const double* cc = P.c.data();
            const double dot = psum(na, [&](size_t j) { return f[j] * (h[j] - cc[j] + al * (u[j] - ld)) * (x[j] - xprev[j]); });
            if (dot > 0) tk = 1.0;
            const double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * tk * tk)), beta = (tk - 1.0) / tn;
            tk = tn;
            pfor(na, [&](size_t a, size_t b) {
                for (size_t j = a; j < b; j++) {
                    const double v = x[j] + beta * (x[j] - xprev[j]);
                    xprev[j] = x[j];
                    u[j] = std::min(P.hi, std::max(P.lo, v));
                }
            });
        } else {
            u.swap(x);
        }
        Q = P.value(u, f, fit);
        const double prevbest = bestQ;
        if (Q < bestQ) { bestQ = Q; if (mode == 1) best = u; }
        if (st.prod_target < 0 && Q <= target) st.prod_target = pc.prod();
        if (rtol > 0 && it >= 2 && prevbest - bestQ < rtol * std::fabs(bestQ)) { it++; break; }
    }
    if (mode == 1 && bestQ < Q) { u = best; Q = bestQ; }
    st.Q = Q; st.iters = it; st.prod = pc.prod();
}

// mmf: mm with momentum as mmacc, the MM point by two Newton steps on h/f e^v + alpha (v - ln m) = c from v = u
// (no logarithm; the Wright omega only where the first step exceeds 0.1 in ln f), the restart sum in the same pass
static void mmf_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st, double rtol) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na;
    std::vector<double> f, fit, h, wf, x(na), om, xprev = u;
    double Q = P.value(u, f, fit);
    st.Q0 = Q;
    double bestQ = Q;
    double tk = 1.0;
    int it = 0;
    const double al = P.alpha, ld = P.ld, la = std::log(P.alpha), lo = P.lo, hi = P.hi;
    const double* cc = P.c.data();
    for (; it < maxit && pc.prod() < budget; it++) {
        hf_of(P, fit, h, wf);
        double dot;
        {
            Timed tm(T_MMMAP);
            if (om.size() != na) om.assign(na, -1.0);
            dot = psum(na, [&](size_t j) {
                const double hj = h[j], uj = u[j], ca = cc[j];
                const double g = hj - ca + al * (uj - ld);   // dQ/df at u
                const double den = hj + al, d1 = -g / den;
                double v;
                if (std::fabs(d1) < 0.1) {
                    const double v1 = uj + d1, e = fexp(d1), he = hj * e;
                    v = v1 - (he + al * (v1 - ld) - ca) / (he + al);
                } else if (hj > 0) {
                    const double t = flog(hj) - uj + ld - la + ca / al;
                    const double o = wright_omega(t, om[j]);
                    om[j] = o;
                    v = ld + ca / al - o;
                } else {
                    v = ld + ca / al;
                }
                const double xj = std::min(hi, std::max(lo, v));
                x[j] = xj;
                return f[j] * g * (xj - xprev[j]);
            });
        }
        if (dot > 0) tk = 1.0;   // restart: the momentum points uphill
        const double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * tk * tk)), beta = (tk - 1.0) / tn;
        tk = tn;
        pfor(na, [&](size_t a, size_t b) {
            for (size_t j = a; j < b; j++) {
                const double v = x[j] + beta * (x[j] - xprev[j]);
                xprev[j] = x[j];
                u[j] = std::min(hi, std::max(lo, v));
            }
        });
        Q = P.value(u, f, fit);
        const double prevbest = bestQ;
        if (Q < bestQ) bestQ = Q;
        if (st.prod_target < 0 && Q <= target) st.prod_target = pc.prod();
        for (int c = 0; c < 6; c++) if (st.Qat[c] == 0 && pc.prod() + 2 > kCheck[c]) st.Qat[c] = Q;
        if (rtol > 0 && it >= 2 && prevbest - bestQ < rtol * std::fabs(bestQ)) { it++; break; }
    }
    st.Q = Q; st.iters = it; st.prod = pc.prod();
}

// the MM point x of u (h = A'(w A f) at u): as in mmf
static inline double mm_point(double hj, double uj, double ca, double al, double ld, double la, double lo, double hi, double& omj) {
    const double g = hj - ca + al * (uj - ld);
    const double den = hj + al, d1 = -g / den;
    double v;
    if (std::fabs(d1) < 0.1) {
        const double v1 = uj + d1, e = fexp(d1), he = hj * e;
        v = v1 - (he + al * (v1 - ld) - ca) / (he + al);
    } else if (hj > 0) {
        const double t = flog(hj) - uj + ld - la + ca / al;
        const double o = wright_omega(t, omj);
        omj = o;
        v = ld + ca / al - o;
    } else {
        v = ld + ca / al;
    }
    return std::min(hi, std::max(lo, v));
}

// mmaa: Anderson acceleration (type II, memory m, MS_SOLVER_AA) of the MM map G in u, least squares in the metric of
// f at the start; a step that raises Q is replaced by the plain MM step and the memory cleared
static void mmaa_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st, double rtol) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na;
    const int m = getenv("MS_SOLVER_AA") ? std::max(1, std::min(10, atoi(getenv("MS_SOLVER_AA")))) : 5;
    std::vector<double> f, fit, h, wf, g(na), om(na, -1.0), Fk(na), Fp(na), gp(na), un, fn, fitn, W;
    std::vector<std::vector<float>> dF, dG;   // ring, oldest first
    std::vector<double> A;                      // <dF_i, dF_j>_W, m x m
    double Q = P.value(u, f, fit);
    W = f;
    st.Q0 = Q;
    const double al = P.alpha, ld = P.ld, la = std::log(P.alpha), lo = P.lo, hi = P.hi;
    const double* cc = P.c.data();
    int it = 0, rejects = 0;
    bool have_prev = false;
    for (; it < maxit && pc.prod() < budget; it++) {
        hf_of(P, fit, h, wf);
        {
            Timed tm(T_MMMAP);
            pfor(na, [&](size_t a, size_t b) {
                for (size_t j = a; j < b; j++) { g[j] = mm_point(h[j], u[j], cc[j], al, ld, la, lo, hi, om[j]); Fk[j] = g[j] - u[j]; }
            });
        }
        if (have_prev) {   // the new differences
            std::vector<float> df(na), dg(na);
            pfor(na, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) { df[j] = (float)(Fk[j] - Fp[j]); dg[j] = (float)(g[j] - gp[j]); } });
            if ((int)dF.size() == m) { dF.erase(dF.begin()); dG.erase(dG.begin()); }
            dF.push_back(std::move(df)); dG.push_back(std::move(dg));
        }
        Fp = Fk; gp = g; have_prev = true;
        const int k = (int)dF.size();
        std::vector<double> gam(k, 0.0);
        if (k > 0) {   // min |Fk - dF gam|_W by the normal equations (all dots each step: k <= m)
            const int nd = k * (k + 1) / 2 + k;
            std::vector<double> d(nd);
            psums(na, nd, d.data(), [&](size_t i0, size_t i1, double* acc) {
                int q = 0;
                for (int a = 0; a < k; a++) {
                    const float* fa = dF[a].data();
                    for (int b = a; b < k; b++) {
                        const float* fb = dF[b].data();
                        double v = 0;
                        for (size_t j = i0; j < i1; j++) v += W[j] * (double)fa[j] * (double)fb[j];
                        acc[q++] = v;
                    }
                    double v = 0;
                    for (size_t j = i0; j < i1; j++) v += W[j] * (double)fa[j] * Fk[j];
                    acc[q++] = v;
                }
            });
            std::vector<double> M((size_t)k * k), rhs(k);
            int q = 0;
            for (int a = 0; a < k; a++) {
                for (int b = a; b < k; b++) M[(size_t)a * k + b] = M[(size_t)b * k + a] = d[q++];
                rhs[a] = d[q++];
            }
            double tr = 0;
            for (int a = 0; a < k; a++) tr += M[(size_t)a * k + a];
            for (int a = 0; a < k; a++) M[(size_t)a * k + a] += 1e-10 * tr + 1e-300;
            // Cholesky
            std::vector<double> L((size_t)k * k, 0.0);
            bool ok = true;
            for (int i = 0; i < k && ok; i++)
                for (int j = 0; j <= i; j++) {
                    double sm = M[(size_t)i * k + j];
                    for (int t = 0; t < j; t++) sm -= L[(size_t)i * k + t] * L[(size_t)j * k + t];
                    if (i == j) { if (!(sm > 0)) { ok = false; break; } L[(size_t)i * k + i] = std::sqrt(sm); }
                    else L[(size_t)i * k + j] = sm / L[(size_t)j * k + j];
                }
            if (ok) {
                std::vector<double> z(k);
                for (int i = 0; i < k; i++) { double sm = rhs[i]; for (int t = 0; t < i; t++) sm -= L[(size_t)i * k + t] * z[t]; z[i] = sm / L[(size_t)i * k + i]; }
                for (int i = k - 1; i >= 0; i--) { double sm = z[i]; for (int t = i + 1; t < k; t++) sm -= L[(size_t)t * k + i] * gam[t]; gam[i] = sm / L[(size_t)i * k + i]; }
                for (double& v : gam) if (!std::isfinite(v)) v = 0;
            } else std::fill(gam.begin(), gam.end(), 0.0);
        }
        un.resize(na);
        pfor(na, [&](size_t a, size_t b) {
            for (size_t j = a; j < b; j++) {
                double v = g[j];
                for (int i = 0; i < k; i++) v -= gam[i] * (double)dG[i][j];
                un[j] = std::min(hi, std::max(lo, v));
            }
        });
        double Qn = P.value(un, fn, fitn);
        if (k > 0 && !(Qn <= Q)) {   // safeguard: the plain MM step (Q falls)
            rejects++;
            dF.clear(); dG.clear(); have_prev = false;
            un = g;
            Qn = P.value(un, fn, fitn);
        }
        u.swap(un); f.swap(fn); fit.swap(fitn);
        const double Qold = Q;
        Q = Qn;
        if (st.prod_target < 0 && Q <= target) st.prod_target = pc.prod();
        if (rtol > 0 && it >= 2 && Qold - Q < rtol * std::fabs(Q)) { it++; break; }
    }
    if (getenv("MS_DEBUG_SOLVER")) fprintf(stderr, "    mmaa: %d steps, %d rejected\n", it, rejects);
    st.Q = Q; st.iters = it; st.prod = pc.prod();
}

// mmlb: mmf for MS_SOLVER_MM1 products (default 20), then lbfgs_b (its usual iterations) from there
// SQUAREM (scheme S3) on the MM map, with a plain step when the extrapolation does not lower Q
static void squarem_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na;
    std::vector<double> f, fit, h, wf, u1, u2, up, om, f1, fit1, fp, fitp;
    double Q = P.value(u, f, fit);
    st.Q0 = Q;
    int it = 0;
    double amax = 1.0;   // the largest step length so far allowed (grows while the extrapolations succeed)
    for (; it < maxit && pc.prod() < budget; it++) {
        hf_of(P, fit, h, wf);
        mm_map(P, u, h, u1, om);
        P.value(u1, f1, fit1);   // (the fit at u1 for its MM step)
        hf_of(P, fit1, h, wf);
        mm_map(P, u1, h, u2, om);
        double rr[2];
        psums(na, 2, rr, [&](size_t i0, size_t i1, double* acc) {
            double a = 0, b = 0;
            for (size_t j = i0; j < i1; j++) { const double r = u1[j] - u[j], v = u2[j] - 2 * u1[j] + u[j]; a += r * r; b += v * v; }
            acc[0] = a; acc[1] = b;
        });
        double al = rr[1] > 0 ? -std::sqrt(rr[0] / rr[1]) : -1.0;
        al = std::max(al, -amax);
        bool ok = false;
        if (al < -1.0) {
            up.resize(na);
            pfor(na, [&](size_t a, size_t b) {
                for (size_t j = a; j < b; j++) {
                    const double r = u1[j] - u[j], v = u2[j] - 2 * u1[j] + u[j];
                    up[j] = std::min(P.hi, std::max(P.lo, u[j] - 2 * al * r + al * al * v));
                }
            });
            const double Qp = P.value(up, fp, fitp);
            if (Qp <= Q) { u.swap(up); f.swap(fp); fit.swap(fitp); Q = Qp; ok = true; amax = std::max(1.0, -al) * 4.0; }
            else amax = std::max(1.0, -al / 2.0);
        }
        if (!ok) { u.swap(u2); Q = P.value(u, f, fit); }
        if (st.prod_target < 0 && Q <= target) st.prod_target = pc.prod();
    }
    st.Q = Q; st.iters = it; st.prod = pc.prod();
}

// truncated Newton in f with the majoriser's diagonal as the preconditioner, steps in u
static void tn_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st, int cgmax) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na, ng = P.ng;
    const double al = P.alpha, ld = P.ld;
    const double* ww = P.w->data(); const double* yy = P.y->data(); const double* cc = P.c.data();
    std::vector<double> f, fit, s(na), g(na), Minv(na), x(na), r(na), z(na), p(na), Ap(na), Apg(ng), un(na), fn, fitn;
    std::vector<char> fix(na);
    double gn0 = 1;
    double Q = P.value(u, f, fit);
    st.Q0 = Q;
    int it = 0;
    for (; it < maxit && pc.prod() < budget; it++) {
        // gradient in f: A'w(fit - y) + alpha (u - ld); s = A'w(fit - y), h = s + c
        {
            std::vector<double> wr(ng);
            pfor(ng, [&](size_t a, size_t b) { for (size_t i = a; i < b; i++) wr[i] = ww[i] * (fit[i] - yy[i]); });
            P.model->adjoint_f(wr.data(), s.data());
        }
        pfor(na, [&](size_t a, size_t b) {
            for (size_t j = a; j < b; j++) {
                g[j] = s[j] + al * (u[j] - ld);
                const double hj = std::max(s[j] + cc[j], 0.0);
                // masses at a bound that the gradient pushes outwards stay there
                fix[j] = (u[j] <= P.lo && g[j] > 0) || (u[j] >= P.hi && g[j] < 0);
                Minv[j] = fix[j] ? 0.0 : f[j] / (hj + al);
            }
        });
        // PCG on (H + alpha / f) x = -g over the free masses
        const double gn = std::sqrt(psum(na, [&](size_t j) { return fix[j] ? 0.0 : g[j] * g[j]; }));
        if (it == 0) gn0 = gn;
        const double eta = std::min(0.5, std::sqrt(gn / std::max(gn0, 1e-300)));
        pfor(na, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) { x[j] = 0; r[j] = fix[j] ? 0.0 : -g[j]; z[j] = Minv[j] * r[j]; p[j] = z[j]; } });
        double rz = psum(na, [&](size_t j) { return r[j] * z[j]; });
        for (int k = 0; k < cgmax; k++) {
            P.model->forward(p.data(), Apg.data());
            pfor(ng, [&](size_t a, size_t b) { for (size_t i = a; i < b; i++) Apg[i] *= ww[i]; });
            P.model->adjoint_f(Apg.data(), Ap.data());
            const double pAp = psum(na, [&](size_t j) {
                if (fix[j]) { Ap[j] = 0; return 0.0; }
                Ap[j] += al * p[j] / f[j];
                return p[j] * Ap[j];
            });
            if (!(pAp > 0)) break;
            const double a_ = rz / pAp;
            double rn2[2];
            psums(na, 2, rn2, [&](size_t i0, size_t i1, double* acc) {
                double s1 = 0, s2 = 0;
                for (size_t j = i0; j < i1; j++) {
                    x[j] += a_ * p[j]; r[j] -= a_ * Ap[j]; z[j] = Minv[j] * r[j];
                    s1 += r[j] * z[j]; s2 += r[j] * r[j];
                }
                acc[0] = s1; acc[1] = s2;
            });
            const double b_ = rn2[0] / rz;
            rz = rn2[0];
            if (std::sqrt(rn2[1]) <= eta * gn) break;
            pfor(na, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) p[j] = z[j] + b_ * p[j]; });
        }
        // the step in u: du = x / f (at most 3 in ln f), backtracking on Q (Armijo)
        std::vector<double> du(na);
        pfor(na, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) du[j] = std::min(3.0, std::max(-3.0, x[j] / f[j])); });
        const double gd = psum(na, [&](size_t j) { return f[j] * g[j] * du[j]; });
        if (!(gd < 0)) break;
        double stp = 1.0, Qn = Q;
        bool ok = false;
        for (int ls = 0; ls < 12; ls++) {
            un.resize(na);
            pfor(na, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) un[j] = std::min(P.hi, std::max(P.lo, u[j] + stp * du[j])); });
            Qn = P.value(un, fn, fitn);
            if (Qn <= Q + 1e-4 * stp * gd) { ok = true; break; }
            stp *= 0.5;
        }
        if (!ok) break;
        u.swap(un); f.swap(fn); fit.swap(fitn); Q = Qn;
        if (st.prod_target < 0 && Q <= target) st.prod_target = pc.prod();
    }
    st.Q = Q; st.iters = it; st.prod = pc.prod();
}

// Newton-CG on the dual: v = W^-1/2 lambda (data space), s = A'(sqrt(w) v), f(v) = clip(m exp(-s / alpha)),
// F(v) = v'v / 2 + b'v - sum (alpha phi(f) + s f), b = sqrt(w) y; grad = v + b - sqrt(w) A f; Hessian
// I + sqrt(w) A diag(f / alpha) A' sqrt(w) (masses inside the box). Start: v = sqrt(w) (A f - y) of the given u.
static void dual_solve(Inner& P, std::vector<double>& u, int maxit, long budget, double target, SolveStats& st, int cgmax) {
    ProdCount pc;
    const size_t na = P.na, ng = P.ng;
    const double al = P.alpha, ld = P.ld, dl = P.dflt;
    const double* ww = P.w->data(); const double* yy = P.y->data();
    std::vector<double> sw(ng), b(ng), v(ng), grad(ng), s(na), lf(na), f(na), fit(ng), tmpg(ng), tmpm(na);
    std::vector<double> x(ng), r(ng), p(ng), Ap(ng), vn(ng), sn(na), lfn(na), fnn(na), fitn(ng), gradn(ng), fdiag(na);
    pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) { sw[i] = std::sqrt(ww[i]); b[i] = sw[i] * yy[i]; } });
    double Qp0;
    {
        std::vector<double> f0, fit0;
        Qp0 = P.value(u, f0, fit0);
        pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) v[i] = sw[i] * (fit0[i] - yy[i]); });
    }
    st.Q0 = Qp0;
    // F, its gradient and the primal image at v (one adjoint, one forward); returns F
    auto evalF = [&](const std::vector<double>& vv, std::vector<double>& ss, std::vector<double>& lff, std::vector<double>& ff,
                     std::vector<double>& ft, std::vector<double>& gr) {
        pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) tmpg[i] = sw[i] * vv[i]; });
        P.model->adjoint_f(tmpg.data(), ss.data());
        const double conj = psum(na, [&](size_t j) {
            const double l = std::min(P.hi, std::max(P.lo, ld - ss[j] / al));
            lff[j] = l;
            const double fh = fexp(l);
            ff[j] = fh;
            return al * (fh * (l - ld) - fh + dl) + ss[j] * fh;
        });
        P.model->forward(ff.data(), ft.data());
        const double q = psum(ng, [&](size_t i) { gr[i] = vv[i] + b[i] - sw[i] * ft[i]; return 0.5 * vv[i] * vv[i] + b[i] * vv[i]; });
        return q - conj;
    };
    auto primalQ = [&](const std::vector<double>& lff, const std::vector<double>& ff, const std::vector<double>& ft) {
        return P.chi2half(ft) + al * P.entropy_term(lff, ff);
    };
    double F = evalF(v, s, lf, f, fit, grad);
    double Qp = primalQ(lf, f, fit);
    double bestQ = Qp;
    std::vector<double> bestu = lf;
    if (Qp0 < bestQ) { bestQ = Qp0; bestu = u; }
    int it = 0;
    double gn0 = 1;
    double lastgap = Qp + F;
    for (; it < maxit && pc.prod() < budget; it++) {
        const double gn = std::sqrt(pdot(grad.data(), grad.data(), ng));
        pfor(na, [&](size_t a, size_t c) { for (size_t j = a; j < c; j++) fdiag[j] = (lf[j] > P.lo && lf[j] < P.hi) ? f[j] / al : 0.0; });
        // CG on H x = -grad
        pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) { x[i] = 0; r[i] = -grad[i]; p[i] = r[i]; } });
        double rr = gn * gn;
        if (it == 0) gn0 = gn;
        const double eta = std::min(0.5, std::sqrt(gn / std::max(gn0, 1e-300)));
        for (int k = 0; k < cgmax; k++) {
            pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) tmpg[i] = sw[i] * p[i]; });
            P.model->adjoint_f(tmpg.data(), tmpm.data());
            pfor(na, [&](size_t a, size_t c) { for (size_t j = a; j < c; j++) tmpm[j] *= fdiag[j]; });
            P.model->forward(tmpm.data(), Ap.data());
            const double pAp = psum(ng, [&](size_t i) { Ap[i] = p[i] + sw[i] * Ap[i]; return p[i] * Ap[i]; });
            if (!(pAp > 0)) break;
            const double a_ = rr / pAp;
            const double rn = psum(ng, [&](size_t i) { x[i] += a_ * p[i]; r[i] -= a_ * Ap[i]; return r[i] * r[i]; });
            const double b_ = rn / rr;
            rr = rn;
            if (std::sqrt(rn) <= eta * gn) break;
            pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) p[i] = r[i] + b_ * p[i]; });
        }
        const double gd = pdot(grad.data(), x.data(), ng);
        if (!(gd < 0)) break;
        double stp = 1.0, Fn = F;
        bool ok = false;
        for (int ls = 0; ls < 20; ls++) {
            pfor(ng, [&](size_t a, size_t c) { for (size_t i = a; i < c; i++) vn[i] = v[i] + stp * x[i]; });
            Fn = evalF(vn, sn, lfn, fnn, fitn, gradn);
            if (std::isfinite(Fn) && Fn <= F + 1e-4 * stp * gd) { ok = true; break; }
            stp *= 0.5;
        }
        if (!ok) break;
        v.swap(vn); s.swap(sn); lf.swap(lfn); f.swap(fnn); fit.swap(fitn); grad.swap(gradn); F = Fn;
        Qp = primalQ(lf, f, fit);
        lastgap = Qp + F;
        if (Qp < bestQ) { bestQ = Qp; bestu = lf; }
        if (st.prod_target < 0 && Qp <= target) st.prod_target = pc.prod();
        if (getenv("MS_DEBUG_SOLVER")) fprintf(stderr, "    dual it %d: F %.8g primal Q %.8g gap %.4g step %.3g prod %ld\n", it, F, Qp, lastgap, stp, pc.prod());
    }
    u = bestu;
    st.Q = bestQ; st.iters = it; st.prod = pc.prod();
    (void)lastgap;
}

// lbmm: L-BFGS (two-loop recursion, MS_SOLVER_M pairs, default 6) in u whose first inverse Hessian is the majoriser's
// curvature at the current point, 1 / (f (H f + alpha)), recomputed every step (times s'y / y'H0 y), masses at a bound
// that the gradient pushes outwards held, the trial points projected onto the box, Armijo backtracking from step 1
static void lbmm_solve(Objective& obj, Inner& P, std::vector<double>& u, long budget, double target, SolveStats& st, double rtol) {
    ProdCount pc;
    P.make_c();
    const size_t na = P.na;
    const int m = getenv("MS_SOLVER_M") ? std::max(1, std::min(20, atoi(getenv("MS_SOLVER_M")))) : 6;
    const double lo = P.lo, hi = P.hi, al = P.alpha;
    const double* cc = P.c.data();
    obj.cvec = &P.c;
    std::vector<double> g, gn, hd(na), d(na), xn(na), q(na);
    std::vector<char> act(na);
    std::vector<std::vector<float>> S, Y;
    std::vector<double> rho;
    double fx = obj(u, g);
    st.Q0 = fx;
    int it = 0;
    double gam = 1.0;
    for (; pc.prod() < budget; it++) {
        // H0 and the active set at u
        pfor(na, [&](size_t a, size_t b) {
            for (size_t j = a; j < b; j++) {
                hd[j] = 1.0 / std::max(obj.f[j] * (std::max(obj.adj[j] + cc[j], 0.0) + al), 1e-300);
                act[j] = (u[j] <= lo && g[j] > 0) || (u[j] >= hi && g[j] < 0);
                q[j] = act[j] ? 0.0 : g[j];
            }
        });
        const int k = (int)S.size();
        std::vector<double> a(k);
        for (int i = k - 1; i >= 0; i--) {
            const float* si = S[i].data(); const float* yi = Y[i].data();
            a[i] = rho[i] * psum(na, [&](size_t j) { return (double)si[j] * q[j]; });
            const double ai = a[i];
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) q[j] -= ai * (double)yi[j]; });
        }
        pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) q[j] *= gam * hd[j]; });
        for (int i = 0; i < k; i++) {
            const float* si = S[i].data(); const float* yi = Y[i].data();
            const double bi = rho[i] * psum(na, [&](size_t j) { return (double)yi[j] * q[j]; });
            const double c = a[i] - bi;
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) q[j] += c * (double)si[j]; });
        }
        // d = -q on the free masses
        double dg = psum(na, [&](size_t j) { d[j] = act[j] ? 0.0 : -q[j]; return d[j] * g[j]; });
        if (!(dg < 0)) {   // not a descent direction: the scaled gradient, memory cleared
            S.clear(); Y.clear(); rho.clear();
            dg = psum(na, [&](size_t j) { d[j] = act[j] ? 0.0 : -hd[j] * g[j]; return d[j] * g[j]; });
            if (!(dg < 0)) break;
        }
        double stp = 1.0, fn = fx;
        bool ok = false;
        for (int ls = 0; ls < 12 && pc.prod() < budget + 40; ls++) {
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) xn[j] = std::min(hi, std::max(lo, u[j] + stp * d[j])); });
            fn = obj(xn, gn);
            const double gd = psum(na, [&](size_t j) { return g[j] * (xn[j] - u[j]); });
            if (std::isfinite(fn) && fn <= fx + 1e-4 * std::min(gd, 0.0) && fn < fx) { ok = true; break; }
            stp *= 0.25;
        }
        if (!ok) {
            if (S.empty()) break;
            S.clear(); Y.clear(); rho.clear();
            continue;
        }
        // the pair
        std::vector<float> sv(na), yv(na);
        double sy[3];
        psums(na, 3, sy, [&](size_t i0, size_t i1, double* acc) {
            double v0 = 0, v1 = 0, v2 = 0;
            for (size_t j = i0; j < i1; j++) {
                const double sj = xn[j] - u[j], yj = gn[j] - g[j];
                sv[j] = (float)sj; yv[j] = (float)yj;
                v0 += sj * yj; v1 += yj * yj * hd[j]; v2 += sj * sj;
            }
            acc[0] = v0; acc[1] = v1; acc[2] = v2;
        });
        if (sy[0] > 1e-12 * std::sqrt(sy[1] * sy[2] + 1e-300) && sy[0] > 0) {
            if ((int)S.size() == m) { S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); }
            S.push_back(std::move(sv)); Y.push_back(std::move(yv)); rho.push_back(1.0 / sy[0]);
            gam = sy[1] > 0 ? sy[0] / sy[1] : 1.0;
        }
        const double fold = fx;
        u.swap(xn); g.swap(gn); fx = fn;
        if (st.prod_target < 0 && fx <= target) st.prod_target = pc.prod();
        for (int c = 0; c < 6; c++) if (st.Qat[c] == 0 && pc.prod() + 2 > kCheck[c]) st.Qat[c] = fx;
        if (rtol > 0 && it >= 2 && fold - fx < rtol * std::fabs(fx)) { it++; break; }
    }
    obj.cvec = nullptr;
    st.Q = fx; st.iters = it; st.prod = pc.prod();
}

// one solve by name (lbfgs, lbfgsmm, mm, mmacc, squarem, tn, dual) at the problem of obj; target: Q to record the
// products for (or -inf); budget: products at most
static SolveStats run_solver(const std::string& name, Objective& obj, std::vector<double>& u, double lo_u, double hi_u,
                             int iterations, long budget, double target, double rtol = 0) {
    SolveStats st;
    Timed tsol(T_SOLVER);
    const auto t0 = std::chrono::steady_clock::now();
    Inner P(obj, lo_u, hi_u);
    const int big = 1000000;
    const int cgmax = getenv("MS_SOLVER_CG") ? atoi(getenv("MS_SOLVER_CG")) : 20;
    if (name == "lbfgs" || name == "lbfgsmm") {
        ProdCount pc;
        if (name == "lbfgsmm") { P.make_c(); obj.cvec = &P.c; }
        lbfgs_b(obj, u, lo_u, hi_u, iterations);
        obj.cvec = nullptr;
        st.prod = pc.prod();
        std::vector<double> f, fit;
        st.Q = P.value(u, f, fit);
    } else if (name == "mm") mm_solve(P, u, big, budget, target, st, 0, rtol);
    else if (name == "mmacc") mm_solve(P, u, big, budget, target, st, 1, rtol);
    else if (name == "mmf") mmf_solve(P, u, big, budget, target, st, rtol);
    else if (name == "mmaa") mmaa_solve(P, u, big, budget, target, st, rtol);
    else if (name == "lbmm") lbmm_solve(obj, P, u, budget, target, st, rtol);
    else if (name == "mmlb") {
        const long b1 = getenv("MS_SOLVER_MM1") ? atol(getenv("MS_SOLVER_MM1")) : 20;
        ProdCount pc;
        mmf_solve(P, u, big, b1, target, st, 0);
        const long p1 = pc.prod();
        lbfgs_b(obj, u, lo_u, hi_u, iterations);
        st.prod = pc.prod();
        (void)p1;
    }
    else if (name == "squarem") squarem_solve(P, u, big, budget, target, st);
    else if (name == "tn") tn_solve(P, u, big, budget, target, st, cgmax);
    else if (name == "dual") dual_solve(P, u, big, budget, target, st, cgmax);
    else throw std::invalid_argument("MS_SOLVER: unknown solver " + name);
    if (getenv("MS_DEBUG") || getenv("MS_SOLVER_CMP")) {   // the gap of the result (not counted)
        const long fsave = prof_n[T_FORWARD], asave = prof_n[T_ADJOINT];
        std::vector<double> f, fit;
        pfor(u.size(), [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) u[j] = std::min(hi_u, std::max(lo_u, u[j])); });
        st.Q = P.value(u, f, fit);
        st.gap = P.gap(u, f, fit);
        st.chi2 = 2.0 * P.chi2half(fit);
        prof_n[T_FORWARD] = fsave; prof_n[T_ADJOINT] = asave;
    }
    st.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return st;
}
