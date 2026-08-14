%% vf_fit_order8.m
% Low-order rational approximation Ga(s) of frequency-response data
% G(j*omg(k)), k=1..M, using (relaxed) Vector Fitting.
%
% Why Vector Fitting and not invfreqs/polynomial LS:
% Your data spans 6 decades (1e-3..1e3 rad/s). Fitting an order-8
% transfer function by working directly with powers of s (as invfreqs
% does) means building s^8 at both ends of the band -- a dynamic range
% of ~1e48 within a single column of the least-squares matrix. This was
% tested directly: it frequently produces unstable or badly wrong poles
% even when the aggregate curve looks fine. Vector Fitting instead uses
% a partial-fraction (pole) basis, which stays well conditioned
% regardless of bandwidth, and is the standard tool in the literature
% for exactly this problem (rational approximation of tabulated/
% infinite-dimensional frequency response data). Reference: Gustavsen &
% Semlyen (1999), "Rational approximation of frequency domain responses
% by vector fitting"; relaxed version: Gustavsen (2006).
%
% -------------------------------------------------------------------
% USAGE
%   omg = logspace(-3,3,2501);
%   G   = ...;                  % complex, G(k) = G(j*omg(k)), same size as omg
%   n     = 8;                  % desired order of Ga(s)
%   niter = 15;                 % vector fitting iterations (15-20 is plenty)
%
%   [num, den, poles, relerr] = vf_fit_order8(omg, G, n, niter);
%   Ga = tf(num, den);          % needs Control System Toolbox (optional)
%
% Everything below is self-contained (no toolbox required for the fit
% itself; only 'tf'/'bode' at the very end are optional convenience calls).
%
% RELATIVE DEGREE: the fitted model is forced to be strictly proper with
% relative degree exactly 1 (numel(den)-numel(num) == 1), i.e. no direct
% feedthrough term. This is enforced structurally (no constant "+d" term
% anywhere in the pole-relocation or residue-fitting steps), not by
% truncating the numerator after the fact, so the pole search itself is
% consistent with the constraint at every iteration.
% -------------------------------------------------------------------

function [num, den, poles, relerr] = vf_fit_order8(omg, G, n, niter)

if nargin < 3 || isempty(n),     n = 8;  end
if nargin < 4 || isempty(niter), niter = 15; end

omg = omg(:);
G   = G(:);
s   = 1i*omg;
wt  = ones(size(s));   % uniform weighting; see note at bottom for alternatives

% ---- initial poles: complex-conjugate pairs, log-spaced across the band
poles = init_poles(omg, n);

% ---- iteratively relocate poles (linear "pole relocating" LS problem)
for it = 1:niter
    poles = relocate_poles_relaxed(s, G, poles, wt);
end

% ---- final residues with poles held fixed (linear LS)
[c, d, ord] = final_residues(s, G, poles, wt);

% ---- convert pole/residue model to polynomial Ga(s) = num(s)/den(s)
[num, den, ~] = realbasis_to_tf(poles, c, d, ord);
num = real(num);
den = real(den);
% d == 0 by construction (strictly proper fit) => the leading (highest
% power of s) coefficient of num is exactly 0 up to roundoff; drop it so
% numel(den) - numel(num) == 1, i.e. relative degree 1 exactly.
num(1) = [];
poles = roots(den);

% ---- diagnostics
Gfit   = polyval(num, s) ./ polyval(den, s);
relerr = norm(abs(Gfit - G)) / norm(abs(G));
reldeg = numel(den) - numel(num);
fprintf('vf_fit_order8: order %d, relative degree = %d (want 1), max real(pole) = %.3g (want <= 0), relative L2 magnitude error = %.4g\n', ...
        n, reldeg, max(real(poles)), relerr);

% ---- validation plot
figure;
subplot(2,1,1)
loglog(omg, abs(G), 'b', omg, abs(Gfit), 'r--', 'LineWidth', 1.2)
xlabel('\omega [rad/s]'); ylabel('|G|'); legend('data','G_a fit','Location','best'); grid on
title(sprintf('Order-%d fit, relative L2 magnitude error = %.3g', n, relerr))

subplot(2,1,2)
semilogx(omg, unwrap(angle(G))*180/pi, 'b', omg, unwrap(angle(Gfit))*180/pi, 'r--', 'LineWidth', 1.2)
xlabel('\omega [rad/s]'); ylabel('phase [deg]'); grid on

end

% =====================================================================
function poles0 = init_poles(omg, n)
    wmin = min(omg); wmax = max(omg);
    npairs = floor(n/2);
    beta = logspace(log10(wmin), log10(wmax), npairs);
    poles0 = [];
    for k = 1:npairs
        b = beta(k);
        a = -1e-2*b;
        poles0 = [poles0, a-1i*b, a+1i*b]; %#ok<AGROW>
    end
    if mod(n,2) == 1
        poles0 = [poles0, -wmax];
    end
end

% =====================================================================
function [Phi, ord] = build_basis(s, poles)
% Real-valued partial-fraction basis for a (possibly complex) pole set.
% Each real pole contributes 1 column; each complex-conjugate pair
% contributes 2 real-valued combination columns (standard VF trick,
% keeps all unknowns real even though the poles are complex).
    poles = poles(:).';
    N = numel(poles);
    used = false(1,N);
    cols = {};
    ord  = {};
    k = 1;
    while k <= N
        if used(k), k = k+1; continue; end
        p = poles(k);
        if abs(imag(p)) < 1e-10*max(1,abs(p))
            cols{end+1} = 1./(s-p);           %#ok<AGROW>
            ord{end+1}  = {'real', p};         %#ok<AGROW>
            used(k) = true;
        else
            pc = conj(p);
            phi_a = 1./(s-p) + 1./(s-pc);
            phi_b = 1i*(1./(s-p) - 1./(s-pc));
            cols{end+1} = phi_a;               %#ok<AGROW>
            cols{end+1} = phi_b;               %#ok<AGROW>
            ord{end+1}  = {'cplx', p, pc};      %#ok<AGROW>
            used(k) = true;
            for j = k+1:N
                if ~used(j) && abs(poles(j)-pc) < 1e-6*max(1,abs(pc))
                    used(j) = true;
                    break
                end
            end
        end
        k = k+1;
    end
    Phi = cat(2, cols{:});
end

% =====================================================================
function [num, den, ps] = realbasis_to_tf(poles, c, d, ord)
% Convert d + sum_k c_k*phi_k(s)  (real-basis form defined by ord) into
% an explicit polynomial ratio num(s)/den(s), den(s) = prod(s-poles).
    res = []; ps = [];
    idx = 1;
    for kk = 1:numel(ord)
        item = ord{kk};
        if strcmp(item{1}, 'real')
            res(end+1) = c(idx);  ps(end+1) = item{2};             %#ok<AGROW>
            idx = idx+1;
        else
            a = c(idx); b = c(idx+1); p = item{2}; pc = item{3};
            res(end+1) = a+1i*b;  ps(end+1) = p;                    %#ok<AGROW>
            res(end+1) = a-1i*b;  ps(end+1) = pc;                   %#ok<AGROW>
            idx = idx+2;
        end
    end
    den = poly(ps);
    num = d*den;
    for kk = 1:numel(ps)
        others = ps([1:kk-1, kk+1:end]);
        term = res(kk)*poly(others);
        if numel(term) < numel(num)
            term = [zeros(1, numel(num)-numel(term)), term];
        end
        num = num + term;
    end
end

% =====================================================================
function poles = relocate_poles_relaxed(s, G, poles, wt)
% One relaxed-Vector-Fitting pole-relocation step (linear LS problem).
% NOTE: f(s) is modeled as strictly proper (sum c_k*phi_k, no "+d" term)
% so that the fitted model has relative degree exactly 1. This is
% enforced here (not just in final_residues) so the pole search is
% consistent with that constraint at every iteration.
    [Phi, ord] = build_basis(s, poles);
    [Ns, N] = size(Phi);
    A_f    = Phi;                          % f(s) = sum c_k*phi_k   (no +d)
    A_sigc = -G(:).*Phi;
    A_sigd = -G(:);
    top = [A_f, A_sigc, A_sigd];           % unknowns: c(N), ctil(N), dtil
    rhs = zeros(Ns,1);
    Aw = top .* wt(:);
    rw = rhs .* wt(:);
    Ari = [real(Aw); imag(Aw)];
    rri = [real(rw); imag(rw)];
    % relaxed non-triviality constraint (Gustavsen 2006), avoids the
    % trivial all-zero solution without forcing sigma's constant term to 1
    extraRow = zeros(1, 2*N+1);
    extraRow(N+1:2*N) = real(sum(Phi,1));
    extraRow(end) = Ns;
    Ari = [Ari; extraRow];
    rri = [rri; Ns];
    sol = Ari \ rri;
    ctil = sol(N+1:2*N);
    dtil = sol(end);

    [numer_poly, ~, ~] = realbasis_to_tf(poles, ctil, dtil, ord);
    numer_poly = real(numer_poly);           % should already be real (real-basis solve); clean roundoff
    new_poles = roots(numer_poly).';
    rp = real(new_poles) > 0;
    new_poles(rp) = -real(new_poles(rp)) + 1i*imag(new_poles(rp));   % reflect unstable poles
    poles = new_poles;
end

% =====================================================================
function [c, d, ord] = final_residues(s, G, poles, wt)
% Final linear LS for residues with poles held fixed.
% Strictly proper by construction: no constant term is fit, so d = 0
% always, which forces numel(den)-numel(num) == 1 (relative degree 1)
% once realbasis_to_tf builds the polynomial pair.
    [Phi, ord] = build_basis(s, poles);
    A  = Phi;
    Aw = A .* wt(:);
    rw = G(:) .* wt(:);
    Ari = [real(Aw); imag(Aw)];
    rri = [real(rw); imag(rw)];
    sol = Ari \ rri;
    c = sol;
    d = 0;
end

% =====================================================================
% NOTES
%  - Weighting: wt=ones(...) gives a good absolute (not necessarily
%    relative) fit, weighting more heavily where |G| is largest. If you
%    instead want good RELATIVE accuracy even where |G| is tiny (e.g.
%    deep in a notch, or far past roll-off), try wt = 1./abs(G) -- but
%    be aware this can hurt conditioning/stability if G decays very
%    fast at the band edges; test both and compare the Bode overlay.
%  - If your G(s) is genuinely delay-dominated (phase keeps winding up
%    roughly linearly and |G| stays O(1) across the whole 6 decades),
%    no order-8 rational function will fit well over the full band --
%    that's a fundamental limitation, not a bug. Consider restricting
%    omg to the band you actually care about, or raising the order.
%  - Stability is enforced at every iteration (RHP poles are reflected
%    to the LHP), so the returned Ga(s) is guaranteed stable.
