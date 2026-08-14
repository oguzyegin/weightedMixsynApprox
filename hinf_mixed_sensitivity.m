function [gamma_opt, omg, C, S, T, per, candidates] = hinf_mixed_sensitivity(...
    nw1, dw1, nw2, dw2, nmd, dmd, mn, no, varargin)
% HINF_MIXED_SENSITIVITY  Toker-Ozbay / Foias-Ozbay-Tannenbaum solver for the
% scalar 2-block H-infinity mixed-sensitivity problem
%
%     gamma_opt = min over stabilizing C of  || [W1*S ; W2*T] ||_inf
%
% for a (possibly infinite-dimensional / delay) plant factored as
%
%     P(s) = Mn(s)*No(s)/Md(s)
%
% where Md is the FINITE Blaschke product for the plant's unstable poles,
% Mn is the remaining inner (all-pass) part (typically carries the delays),
% and No is the outer (stable, minimum-phase) part.
%
% This consolidates and corrects the logic in the original HINFCON package
% (Toker & Ozbay, 1996) into a single callable function: no more manual
% gamma bisection, no more manually re-running opt.m/opt2.m by hand.
%
% -------------------------------------------------------------------
% INPUTS
%   nw1,dw1   : W1(s) = polyval(nw1,s)/polyval(dw1,s)   (performance weight on S)
%   nw2,dw2   : W2(s) = polyval(nw2,s)/polyval(dw2,s)   (performance weight on T)
%   nmd,dmd   : Md(s) = polyval(nmd,s)/polyval(dmd,s)   (finite Blaschke product;
%               alpha = roots(nmd) are taken as the plant's unstable poles)
%   mn        : function handle, Mn(s)  (inner part not captured by Md; may
%               contain delays, e.g. mn = @(s) exp(-0.2*s).*(...)./(...))
%   no        : function handle, No(s)  (outer, stable min-phase part; may
%               also contain delays)
%
% OPTIONAL NAME-VALUE ARGUMENTS
%   'gmin', 'gmax'   : search bracket for gamma (default 1e-3, 20)
%   'Ncoarse'        : number of coarse grid points for the initial gamma
%                      scan (default 4000)
%   'lw1','lw2'      : log10(freq) range for the returned frequency response
%                      (default -4, 4)
%   'Lpoints'        : number of frequency points (default 2500)
%   'EPS'            : tolerance for classifying roots as RHP/on-axis (default 1e-6)
%
% OUTPUTS
%   gamma_opt  : optimal H-infinity level (the LARGEST gamma at which the
%                Pick matrix is singular -- see Toker & Ozbay 1995/1996)
%   omg        : frequency grid (rad/s), log-spaced
%   C          : optimal controller frequency response, C(1j*omg)
%   S, T       : closed-loop sensitivity/complementary sensitivity, 1/(1+P*C)
%                and P*C/(1+P*C)
%   per        : sqrt(|W1.*S|^2+|W2.*T|^2) -- should be FLAT and equal to
%                gamma_opt at every frequency for a correctly-found optimum
%                (the "all-pass completion" signature of true optimality --
%                always sanity-check this before trusting gamma_opt!)
%   candidates : all near-singular gamma values found during the scan (for
%                diagnostics -- gamma_opt is max(candidates))
%
% -------------------------------------------------------------------
% EXAMPLE (the plant/weights this was validated against):
%
%   nw1 = 0.5; dw1 = [1, 0.00001];
%   nw2 = [0.01, 0.3, 0.2]; dw2 = 1;
%   nmd = [1, -0.468143990182247, 1.233838206204343];
%   dmd = [1,  0.468143990182247, 1.233838206204343];
%   mn  = @(s) exp(-0.2*s).*(5-s.*exp(-0.1*s))./(s+5*exp(-0.1*s));
%   no  = @(s) ((s.^2-0.468143990182247*s+1.233838206204343)./ ...
%               (s.^2+0.468143990182247*s+1.233838206204343)) .* ...
%              ((s+5*exp(-0.1*s))./(3*s+0.5)) ./ ...
%              (1+((2*s+7)./((s+1).*(3*s+0.5))).*exp(-1.5*s) + ...
%                 ((s-1)./((s+1).*(3*s+0.5))).*exp(-2*s)) ./ (s+1);
%
%   [gamma_opt, omg, C, S, T, per] = hinf_mixed_sensitivity(nw1,dw1,nw2,dw2,nmd,dmd,mn,no);
%   % gamma_opt should come out to 1.0025928, per should be flat at that value.
%
% -------------------------------------------------------------------
% VALIDATION: this function has been checked against SEVEN independent
% results and matches every one of them (see test_hinf_mixed_sensitivity.m):
%   - the real HINFCON/MATLAB run above (gamma_opt = 1.0025928)
%   - six worked examples from Ozbay, Gumussoy, Kashima, Yamamoto,
%     "Frequency Domain Techniques for H-infinity Control of Distributed
%     Parameter Systems", SIAM, 2018, Sections 6.3.1-6.3.3:
%       Sec 6.3.1, unstable retarded system   -> gamma_opt = 17.846
%       Sec 6.3.1, AQM/TCP flow control        -> gamma_opt = 171.33
%       Sec 6.3.2, infinitely many zeros in C+ -> gamma_opt = 0.845663
%       Sec 6.3.2, infinitely many poles in C+ -> gamma_opt = 0.9107  (via duality)
%       Sec 6.3.3, fractional-order (weights 1)-> gamma_opt = 1.47023
%       Sec 6.3.3, fractional-order (weights 2)-> gamma_opt = 2.5203
% -------------------------------------------------------------------
% IMPLEMENTATION NOTE ON A BUG THIS FIXES:
% The original HINFCON spec.m is declared as
%     function [nump,denp] = spectral_factorization(nums,dens)
% and is called in opt.m/opt2.m as  [dG,nG] = spec(A,B)
% i.e. dG receives the SCALED factor of the FIRST argument (A), and nG
% receives the PLAIN factor of the SECOND argument (B). A transposed/swapped
% return here silently produces a Pick matrix that is never singular at the
% true optimum -- this was diagnosed and fixed in specfact() below.
% -------------------------------------------------------------------

p = inputParser;
addParameter(p, 'gmin', 1e-3);
addParameter(p, 'gmax', 20.0);
addParameter(p, 'Ncoarse', 4000);
addParameter(p, 'lw1', -4);
addParameter(p, 'lw2', 4);
addParameter(p, 'Lpoints', 2500);
addParameter(p, 'EPS', 1e-6);
parse(p, varargin{:});
opts = p.Results;

nw1 = nw1(:).'; dw1 = dw1(:).';
nw2 = nw2(:).'; dw2 = dw2(:).';
nmd = nmd(:).'; dmd = dmd(:).';

if ~isa(mn, 'function_handle')
    error('hinf_mixed_sensitivity:badMn', ...
        ['mn must be a function handle, e.g. mn = @(s) exp(-0.2*s).*(...)./(...);  ' ...
         'got class %s instead. (If you evaluated it at specific s-values and passed ' ...
         'the resulting numeric array, that is the bug -- pass the @(s)... handle itself.)'], ...
        class(mn));
end
if ~isa(no, 'function_handle')
    error('hinf_mixed_sensitivity:badNo', ...
        'no must be a function handle, e.g. no = @(s) (...); got class %s instead.', class(no));
end

alpha = roots(nmd);   % plant's unstable (RHP) poles, from Md's zeros
md = @(s) polyval(nmd, s)./polyval(dmd, s);

% ---------------- 1. coarse scan for candidate singular gammas ----------------
gam = linspace(opts.gmin, opts.gmax, opts.Ncoarse);
vals = nan(size(gam));
for i = 1:numel(gam)
    v = smin_at(gam(i));
    if ~isempty(v)
        vals(i) = v;
    end
end

candidates = [];
for i = 2:numel(gam)-1
    if isnan(vals(i)); continue; end
    L = vals(i-1); if isnan(L); L = Inf; end
    R = vals(i+1); if isnan(R); R = Inf; end
    if vals(i) <= L && vals(i) <= R && vals(i) < 1e-2
        candidates(end+1) = gam(i); %#ok<AGROW>
    end
end

if isempty(candidates)
    error('hinf_mixed_sensitivity:noSolution', ...
        'No near-singular gamma found in [%.4g, %.4g]. Widen gmin/gmax.', opts.gmin, opts.gmax);
end

% ---------------- 2. refine each candidate (golden-section-style bisection) ----------------
% NOTE: uses candIdx (not k) as the loop variable -- MATLAB/Octave nested
% functions share the ENTIRE parent workspace, including loop variables, so
% a generic name like "k" here would silently collide with any nested
% function that also uses "k" internally (which build_pick_matrix does,
% for its own bookkeeping loops), corrupting this loop's counter mid-flight.
refined = zeros(size(candidates));
for candIdx = 1:numel(candidates)
    g0 = candidates(candIdx);
    lo = max(opts.gmin, g0*0.98);
    hi = min(opts.gmax, g0*1.02);
    for it = 1:60
        m1 = lo + (hi-lo)/3;
        m2 = hi - (hi-lo)/3;
        v1 = smin_at(m1); if isempty(v1); v1 = Inf; end
        v2 = smin_at(m2); if isempty(v2); v2 = Inf; end
        if v1 < v2
            hi = m2;
        else
            lo = m1;
        end
        if hi - lo < 1e-10
            break;
        end
    end
    refined(candIdx) = (lo+hi)/2;
end

% ---------------- 3. gamma_opt = LARGEST singular gamma (Toker-Ozbay) ----------------
gamma_opt = max(refined);
candidates = sort(refined);

% ---------------- 4. build the optimal controller at gamma_opt ----------------
% NOTE: the refined gamma_opt sits exactly on the boundary where the Pick
% matrix becomes singular -- at that exact floating-point value, root-finding
% can occasionally flip to a "non-generic" root count (especially when there
% are few/no unstable poles, e.g. l=0). If that happens, nudge gamma by a
% tiny amount in either direction until a generic point is found, and use
% that (report the nudge to the user so it's never silent).
[M, nF, dF, N] = build_pick_matrix(gamma_opt);
if isempty(M)
    nudged = false;
    for delta = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
        for sgn = [1, -1]
            g_try = gamma_opt + sgn*delta*max(1, gamma_opt);
            [M_try, nF_try, dF_try, N_try] = build_pick_matrix(g_try);
            if ~isempty(M_try)
                fprintf(['hinf_mixed_sensitivity: gamma_opt=%.10g was exactly on a ' ...
                    'non-generic boundary (common when the plant has few/no unstable ' ...
                    'poles); nudged to %.10g to build the controller.\n'], gamma_opt, g_try);
                gamma_opt = g_try;
                M = M_try; nF = nF_try; dF = dF_try; N = N_try;
                nudged = true;
                break
            end
        end
        if nudged; break; end
    end
    if isempty(M)
        error('hinf_mixed_sensitivity:nonGenericAtOptimum', ...
            ['Pick matrix is non-generic at and near gamma_opt=%.10g and could not be ' ...
             'nudged to a generic point. Try a finer Ncoarse, or a tighter [gmin,gmax] ' ...
             'bracket around this value.'], gamma_opt);
    end
end
a1 = polyadd(conv(nw1, starpoly(nw1)), -gamma_opt^2*conv(dw1, starpoly(dw1)));
b1 = gamma_opt^2*conv(dw1, starpoly(dw1));
nE = a1; dE = b1;

[Muu, Suu, Vuu] = svd(M); %#ok<ASGLU>
LL = Vuu(:, size(M,1));
nL = real(LL(N+1:2*N));
dL = real(LL(1:N));

lws = linspace(opts.lw1, opts.lw2, opts.Lpoints+1);
omg = 10.^lws;
s = 1i*omg;
c1 = polyval(nE, s)./polyval(dE, s);
c2 = polyval(nF, s)./polyval(dF, s);
c3 = polyval(nL, s)./polyval(dL, s);
C = c1.*c2.*c3./(1+mn(s).*c2.*c3).*md(s)./no(s);

P = mn(s).*no(s)./md(s);
S = 1./(1+P.*C);
T = 1 - S;
W1v = polyval(nw1, s)./polyval(dw1, s);
W2v = polyval(nw2, s)./polyval(dw2, s);
per = sqrt(abs(W1v.*S).^2 + abs(W2v.*T).^2);

% sanity warning if the "optimal" solution isn't actually flat
if (max(per)-min(per)) > 1e-3*gamma_opt
    warning('hinf_mixed_sensitivity:notFlat', ...
        ['per is not flat (min=%.6g, max=%.6g) -- this gamma_opt may not be ' ...
         'the true optimum. Check gmin/gmax bracket and EPS.'], min(per), max(per));
end

% ===================================================================
% ---- nested helper functions (share nw1,dw1,nw2,dw2,alpha,mn via closure) ----

    function v_ = smin_at(gamma)
        M_ = build_pick_matrix(gamma);
        if isempty(M_)
            v_ = [];
            return
        end
        sv_ = svd(M_);
        v_ = sv_(end);
    end

    function [nF_, dF_, ok] = pick_matrix_parts(gamma)
        a1_ = polyadd(conv(nw1, starpoly(nw1)), -gamma^2*conv(dw1, starpoly(dw1)));
        b1_ = gamma^2*conv(dw1, starpoly(dw1));
        a2_ = polyadd(conv(nw2, starpoly(nw2)), -gamma^2*conv(dw2, starpoly(dw2)));
        b2_ = gamma^2*conv(dw2, starpoly(dw2));
        A1_ = conv(a1_, a2_); B1_ = conv(b1_, b2_);
        A_ = polyadd(B1_, -A1_); B_ = B1_;
        try
            [dG_, nG_] = specfact(A_, B_);
        catch
            nF_ = []; dF_ = []; ok = false;
            return
        end
        eta_ = roots(dw1);
        nF_ = conv(nG_, poly(-eta_));
        dF_ = conv(dG_, poly(eta_));
        ok = true;
    end

    function [M_, nF_, dF_, N_] = build_pick_matrix(gamma)
        [nF_, dF_, ok] = pick_matrix_parts(gamma);
        if ~ok
            M_ = []; N_ = 0;
            return
        end
        a1_ = polyadd(conv(nw1, starpoly(nw1)), -gamma^2*conv(dw1, starpoly(dw1)));
        eta_ = roots(dw1);
        n1_ = numel(eta_);
        beta_ = roots(a1_);
        beta_ = beta_((real(beta_) > opts.EPS) | ...
                      ((abs(real(beta_)) <= opts.EPS) & (imag(beta_) > 0)));
        l_ = numel(alpha);
        if numel(beta_) ~= n1_
            M_ = []; N_ = 0;
            return
        end
        N_ = n1_ + l_;
        rows_ = zeros(2*N_, 2*N_);
        r_ = 1;
        for k_ = 1:numel(beta_)
            pt_ = beta_(k_);
            powers_ = pt_.^(N_-1:-1:0);
            val_ = mn(pt_)*polyval(nF_,pt_)/polyval(dF_,pt_);
            rows_(r_,:) = [powers_, val_*powers_]; r_ = r_+1;
        end
        for k_ = 1:numel(alpha)
            pt_ = alpha(k_);
            powers_ = pt_.^(N_-1:-1:0);
            val_ = mn(pt_)*polyval(nF_,pt_)/polyval(dF_,pt_);
            rows_(r_,:) = [powers_, val_*powers_]; r_ = r_+1;
        end
        for k_ = 1:numel(beta_)
            pt_ = beta_(k_);
            powers_ = (-pt_).^(N_-1:-1:0);
            val_ = mn(pt_)*polyval(nF_,pt_)/polyval(dF_,pt_);
            rows_(r_,:) = [val_*powers_, powers_]; r_ = r_+1;
        end
        for k_ = 1:numel(alpha)
            pt_ = alpha(k_);
            powers_ = (-pt_).^(N_-1:-1:0);
            val_ = mn(pt_)*polyval(nF_,pt_)/polyval(dF_,pt_);
            rows_(r_,:) = [val_*powers_, powers_]; r_ = r_+1;
        end
        M_ = rows_;
    end

end % main function


% =====================================================================
function q = starpoly(p)
% q(s) = p(-s), for a real-coefficient polynomial p given in descending
% power order (MATLAB convention).
n = numel(p);
k = 1:n;
q = p .* (-1).^(n+k);
end

% =====================================================================
function s = polyadd(p, q)
% Add two polynomials (descending power order), right-aligned (i.e.
% aligned at the constant term), trimming leading (high-order) zeros.
EPS = 1e-20;
np = numel(p); nq = numel(q);
n = max(np, nq);
pp = [zeros(1, n-np), p];
qq = [zeros(1, n-nq), q];
s = pp + qq;
while numel(s) > 1 && abs(s(1)) < EPS
    s = s(2:end);
end
end

% =====================================================================
function q = pspecfun(p)
% Build the polynomial from the strictly-left-half-plane roots of p, plus
% one representative from each imaginary-axis conjugate pair.
EPS = 1e-6;
z = roots(p);
zn = z(real(z) < -EPS);
zi = z(abs(real(z)) <= EPS);
zx = [];
while ~isempty(zi)
    zx = [zx; zi(1)]; %#ok<AGROW>
    lz = numel(zi);
    if lz > 1
        [~, n] = min(abs(zi(2:lz) - zi(1)*ones(lz-1,1)));
        n = n + 1;
        zi = zi([2:n-1, n+1:lz]);
    else
        zi = [];
    end
end
q = poly([zn; zx]);
end

% =====================================================================
function [dG, nG] = specfact(nums, dens)
% Spectral factorization of a nonnegative para-Hermitian rational function
% A(s)/B(s) = nums(s)/dens(s). Returns (dG,nG) such that:
%   dG is the SCALED spectral factor of nums (the "A" argument)
%   nG is the PLAIN spectral factor of dens (the "B" argument)
% matching the calling convention [dG,nG] = specfact(A,B) used above.
%
% NOTE: this fixes a return-order bug present in the original HINFCON
% spec.m when read together with how opt.m/opt2.m call it -- see the
% header comment of hinf_mixed_sensitivity.m for the full explanation.
numK = nums(1)/dens(1) * (-1)^((numel(nums)-numel(dens))/2);
if numK < 0
    error('specfact:invalid', 'Cannot do spectral factorization (numK<0).');
end
numK = sqrt(numK);
dG = numK * pspecfun(nums);
nG = pspecfun(dens);
dG = real(dG);
nG = real(nG);
end
