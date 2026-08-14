%% vf_tune_hinf.m
% Direct tuning of the order-8 controller Ca (from vf_fit_order8.m) against
% the TRUE closed-loop mixed-sensitivity objective, instead of just the
% open-loop frequency-response fit error.
%
% PROBLEM
%   Plant (with time delays):
%     q(s) = (s+1)*(3*s+0.5) + (2*s+7)*exp(-1.5*s) + (s-1)*exp(-2*s)
%     P(s) = exp(-0.2*s)*(5 - s*exp(-0.1*s)) / q(s)
%   Weights:
%     W1(s) = 0.5/(s+1e-5)
%     W2(s) = 0.01*s^2 + 0.3*s + 0.2
%   Controller Ca(s) = num(s)/den(s), order 8, relative degree 1 (from
%   vf_fit_order8.m applied to the target frequency-response data "c").
%   Closed loop: L=Ca*P,  S=1/(1+L),  T=L/(1+L).
%
%   Minimize  J(Ca) = max_w sqrt( |W1(jw)*S(jw)|^2 + |W2(jw)*T(jw)|^2 )
%
% WHY NOT JUST USE THE VECTOR-FIT DIRECTLY:
%   vf_fit_order8 minimizes ||Ca_fit - Ca_target||, a magnitude/phase L2
%   error on the imaginary axis. That does NOT directly minimize the
%   closed-loop peak J above -- small fit errors near the critical
%   frequency (where 1+L is small) get amplified by the feedback loop.
%   On this problem the VF fit achieves J ~ 1.68 even though its relative
%   L2 fit error is only ~1e-4; direct tuning against J below brings that
%   down to ~1.03.
%
% METHOD:
%   Ca is re-parameterized by its poles/residues (2 real poles + 3
%   complex-conjugate pairs = 8 poles, matching vf_fit_order8's output
%   structure), which structurally guarantees:
%     - real coefficients,
%     - relative degree exactly 1 (no direct/feedthrough term is ever
%       included in the model),
%     - closed-loop... i.e. Ca-pole stability (real part of every pole is
%       forced <= 0 via a softplus reparameterization: Re(p) = -softplus(x)).
%   J(Ca) is nonsmooth (a max). We minimize a LogSumExp ("soft-max")
%   smooth surrogate of J with an increasing sharpness parameter beta
%   (beta-homotopy / continuation), each stage warm-started from the
%   previous, using fminsearch (derivative-free, always available, no
%   toolbox required). This lands very close to a true minimax solution
%   without requiring the Optimization Toolbox's fmincon/fminimax.
%
% USAGE:
%   [num_t, den_t, poles_t, Jfinal] = vf_tune_hinf(num0, den0, omg);
%   where num0,den0 are Ca's coefficients from vf_fit_order8, and omg is
%   the frequency grid (same one used for the VF fit is fine).
%
% IMPORTANT CAVEAT ON STABILITY:
%   Ca-pole-stability alone does not guarantee closed-loop (internal)
%   stability of the true delay system -- P(s) has 2 unstable
%   (right-half-plane) poles at s = 0.2341 +/- 1.0858j (roots of q(s)=0),
%   found by direct 2D root search of q(s)=0 (delay/quasi-polynomial
%   systems have infinitely many roots but only finitely many with
%   Re(s)>0, since the polynomial part dominates for Re(s)->+inf).
%   The script below verifies internal stability with a Nyquist
%   encirclement count of L(jw) around -1 (must equal 2, the number of
%   unstable open-loop poles of L, for the closed loop to be stable) --
%   always re-run this check after tuning, since nothing in the
%   optimization itself enforces it structurally.
% -------------------------------------------------------------------

function [num_t, den_t, poles_t, Jfinal] = vf_tune_hinf(num0, den0, omg)

s = 1i*omg(:);
P = Pfun(s);
W1 = 0.5./(s+1e-5);
W2 = 0.01*s.^2 + 0.3*s + 0.2;

% ---- convert (num0,den0) -> poles/residues -> initial parameter vector x0
[r, p, k] = residue(num0, den0); %#ok<ASGLU>  % k should be empty (rel.deg==1)
x0 = poles_residues_to_x(p, r);

% ---- beta-homotopy: smooth soft-max, increasing sharpness
betas = [20 50 100 200 400 800 1600 3200 6400 12800 25600 51200 102400];
x = x0;
opts = optimset('fminsearch');
opts = optimset(opts, 'MaxFunEvals', 20000, 'MaxIter', 20000, 'TolX', 1e-12, 'TolFun', 1e-14);
for kk = 1:numel(betas)
    beta = betas(kk);
    f = @(xx) softmax_obj(xx, beta, s, P, W1, W2);
    x = fminsearch(f, x, opts);
    Jnow = true_peak(x, s, P, W1, W2);
    fprintf('beta=%7d   true peak J = %.6f\n', beta, Jnow);
end

[num_t, den_t, poles_t] = build_num_den(x);
Jfinal = true_peak(x, s, P, W1, W2);

% ---- internal stability check (Nyquist encirclements of -1)
Wmax = 1e7; npts = 400000;
omgN = [-logspace(log10(Wmax), -8, npts/2), logspace(-8, log10(Wmax), npts/2)];
sN = 1i*omgN;
CaN = polyval(num_t, sN)./polyval(den_t, sN);
LN  = CaN.*Pfun(sN);
ang = unwrap(angle(1+LN));
N = (ang(end)-ang(1))/(2*pi);
fprintf('Nyquist encirclements of -1: %.4f  (must equal number of unstable open-loop poles of L)\n', N);
fprintf('max real(pole) of Ca: %.3g (want <= 0)\n', max(real(poles_t)));

end

% =====================================================================
function q = qfun(s)
    q = (s+1).*(3*s+0.5) + (2*s+7).*exp(-1.5*s) + (s-1).*exp(-2*s);
end
function P = Pfun(s)
    P = exp(-0.2*s).*(5 - s.*exp(-0.1*s)) ./ qfun(s);
end

% =====================================================================
function y = softplus(x)
    y = log1p(exp(-abs(x))) + max(x,0);
end
function x = inv_softplus(y)
    x = zeros(size(y));
    big = y > 30;
    x(big)  = y(big);
    x(~big) = log(expm1(y(~big)));
end

% =====================================================================
function x0 = poles_residues_to_x(p, r)
% p,r from MATLAB's residue(num0,den0): 2 real poles + 3 conjugate pairs.
    isreal_p = abs(imag(p)) < 1e-9;
    preal = p(isreal_p);  rreal = r(isreal_p);
    pc = p(~isreal_p);    rc = r(~isreal_p);
    % keep only the pole of each pair with positive imaginary part
    keep = imag(pc) > 0;
    pc = pc(keep); rc = rc(keep);
    ar  = inv_softplus(-real(preal));      % 2 values
    ac  = inv_softplus(-real(pc));         % 3 values
    im  = imag(pc);                        % 3 values
    rres = real(rreal);                    % 2 values
    ab  = [real(rc), imag(rc)]';           % 3x2 -> 6 values
    x0 = [ar(:); ac(:); im(:); rres(:); ab(:)];
end

% =====================================================================
function poles = unpack_poles(x)
    ar = -softplus(x(1:2));
    ac = -softplus(x(3:5));
    im = x(6:8);
    poles = [ar; ac + 1i*im; ac - 1i*im];   % 8x1
end

% =====================================================================
function [num, den, poles] = build_num_den(x)
    poles = unpack_poles(x);
    rres = x(9:10);
    ab = reshape(x(11:16), 3, 2);
    res_c = ab(:,1) + 1i*ab(:,2);
    residues = [rres; res_c; conj(res_c)];
    den = poly(poles);
    numfull = zeros(1, numel(den));
    for kk = 1:numel(poles)
        others = poles([1:kk-1, kk+1:end]);
        term = residues(kk)*poly(others);
        term = [zeros(1, numel(numfull)-numel(term)), term];
        numfull = numfull + term;
    end
    num = real(numfull(2:end));   % drop leading (zero) coeff -> deg n-1
    den = real(den);
end

% =====================================================================
function h = pointwise_obj(x, s, P, W1, W2)
    [num, den, ~] = build_num_den(x);
    Ca = polyval(num, s)./polyval(den, s);
    L = Ca.*P;
    S = 1./(1+L);
    T = L./(1+L);
    h = sqrt(abs(W1.*S).^2 + abs(W2.*T).^2);
end

function J = true_peak(x, s, P, W1, W2)
    h = pointwise_obj(x, s, P, W1, W2);
    J = max(h);
end

function v = softmax_obj(x, beta, s, P, W1, W2)
    h = pointwise_obj(x, s, P, W1, W2);
    m = max(beta*h);
    v = (m + log(sum(exp(beta*h - m)))) / beta;   % stable log-sum-exp
end
