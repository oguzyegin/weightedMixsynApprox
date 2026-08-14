function [N,D,tau_n,tau_d,T0,info] = nPdP2YALTA(nP,dP,tol)
%NPD P2YALTA
%
% Convert
%
%                 nP(s)
%        P(s) = ----------
%                 dP(s)
%
% into
%
%                       Nqp(s)
%   P(s) = exp(-T0*s) ----------
%                       Dqp(s)
%
% where
%
%   T0 = nP.OutputDelay - dP.OutputDelay
%
% and
%
%   Nqp(s) = sum_k Nk(s) exp(-tau_n(k)*s)
%   Dqp(s) = sum_k Dk(s) exp(-tau_d(k)*s)
%
% The OutputDelay contribution is extracted explicitly and is NOT
% included in tau_n or tau_d.
%
% Numerator and denominator may have different fundamental internal
% delays.
%
%
% OUTPUTS
% -------
% N
%   Numerator polynomial coefficient matrix.
%   Row k corresponds to tau_n(k).
%
% D
%   Denominator polynomial coefficient matrix.
%   Row k corresponds to tau_d(k).
%
% tau_n
%   Numerator internal-delay grid.
%
% tau_d
%   Denominator internal-delay grid.
%
% T0
%   Overall pure delay:
%
%       T0 = nP.OutputDelay - dP.OutputDelay
%
% info
%   Diagnostic information.
%
%
% Requires:
%   Control System Toolbox
%   Symbolic Math Toolbox
%
% ---------------------------------------------------------------

    if nargin < 3 || isempty(tol)
        tol = 1e-9;
    end

    % ============================================================
    % Basic checks
    % ============================================================

    if ~issiso(nP) || ~issiso(dP)
        error('nP and dP must both be SISO.');
    end

    if nP.Ts ~= 0 || dP.Ts ~= 0
        error('Only continuous-time systems are supported.');
    end


    % ============================================================
    % Extract common/external delay DIRECTLY from OutputDelay
    % ============================================================

    Tn_out = double(nP.OutputDelay);
    Td_out = double(dP.OutputDelay);

    if numel(Tn_out) ~= 1 || numel(Td_out) ~= 1
        error('nP.OutputDelay and dP.OutputDelay must be scalar.');
    end

    T0 = Tn_out - Td_out;

    if abs(T0) < tol
        T0 = 0;
    end


    % ============================================================
    % Remove OutputDelay before reconstructing internal QPs
    % ============================================================

    nP0 = nP;
    dP0 = dP;

    nP0.OutputDelay = 0;
    dP0.OutputDelay = 0;


    % ============================================================
    % Convert each object to rational quasipolynomial form
    %
    %   nP0 = NN / DN
    %   dP0 = ND / DD
    %
    % OutputDelay is NOT included here.
    % ============================================================

    RN = systemToRawQP(nP0,tol);
    RD = systemToRawQP(dP0,tol);


    % ============================================================
    % Form
    %
    %              NN/DD
    %       P0 = ----------
    %              DN/ND
    %
    %          = NN*DD / (DN*ND)
    %
    % with
    %
    %       P = exp(-T0*s) * P0
    % ============================================================

    [numPoly,numDelay] = multiplyQP( ...
        RN.numPoly,RN.numDelay, ...
        RD.denPoly,RD.denDelay,tol);

    [denPoly,denDelay] = multiplyQP( ...
        RN.denPoly,RN.denDelay, ...
        RD.numPoly,RD.numDelay,tol);


    % ============================================================
    % Combine equal-delay terms
    % ============================================================

    [numPoly,numDelay] = combineEqualDelays( ...
        numPoly,numDelay,tol);

    [denPoly,denDelay] = combineEqualDelays( ...
        denPoly,denDelay,tol);


    % ============================================================
    % Remove zero terms
    % ============================================================

    [numPoly,numDelay] = removeZeroTerms( ...
        numPoly,numDelay,tol);

    [denPoly,denDelay] = removeZeroTerms( ...
        denPoly,denDelay,tol);

    if isempty(numPoly)
        error('Resulting numerator is identically zero.');
    end

    if isempty(denPoly)
        error('Resulting denominator is identically zero.');
    end


    % ============================================================
    % Internal quasipolynomials should normally start at delay 0
    %
    % We do NOT use these minimum delays to modify T0.
    % T0 came only from OutputDelay.
    % ============================================================

    minNumDelay = min(numDelay);
    minDenDelay = min(denDelay);

    if minNumDelay > 100*tol
        warning(['Numerator internal quasipolynomial has no ', ...
                 'zero-delay term. Minimum internal delay = %.12g'], ...
                 minNumDelay);
    end

    if minDenDelay > 100*tol
        warning(['Denominator internal quasipolynomial has no ', ...
                 'zero-delay term. Minimum internal delay = %.12g'], ...
                 minDenDelay);
    end

    numDelay(abs(numDelay)<tol) = 0;
    denDelay(abs(denDelay)<tol) = 0;


    % ============================================================
    % Fundamental delays independently
    % ============================================================

    hn = delayGCD(numDelay,tol);
    hd = delayGCD(denDelay,tol);


    % ============================================================
    % Put numerator and denominator on their own regular grids
    % ============================================================

    [N,tau_n,kN] = qpToGrid( ...
        numPoly,numDelay,hn,tol);

    [D,tau_d,kD] = qpToGrid( ...
        denPoly,denDelay,hd,tol);


    % ============================================================
    % Normalize denominator zero-delay polynomial
    % ============================================================

    scale = firstNonzeroCoefficient(D(1,:),tol);

    if ~isempty(scale)
        N = N/scale;
        D = D/scale;
    end

    N(abs(N)<tol) = 0;
    D(abs(D)<tol) = 0;

    % Remove zero numerator rows
    keepN = any(abs(N) > tol,2);
    N = N(keepN,:);
    tau_n = tau_n(keepN);
    
    % Remove zero denominator rows
    keepD = any(abs(D) > tol,2);
    D = D(keepD,:);
    tau_d = tau_d(keepD);


    % ============================================================
    % Diagnostics
    % ============================================================

    info = struct;

    info.TnOutput = Tn_out;
    info.TdOutput = Td_out;
    info.T0       = T0;

    info.hn = hn;
    info.hd = hd;

    info.kN = kN;
    info.kD = kD;

    info.rawNumeratorDelays   = numDelay;
    info.rawDenominatorDelays = denDelay;

    info.nP = RN;
    info.dP = RD;

end


% =================================================================
% SYSTEM -> RAW RATIONAL QUASIPOLYNOMIAL
% =================================================================
function R = systemToRawQP(sys,tol)
%
% Returns
%
%          N(s)
%   sys = ------
%          D(s)
%
% where numerator and denominator are each represented as
%
%       sum_i p_i(s) exp(-tau_i*s)
%
% OutputDelay must already have been removed by nPdP2YALTA.
%

    if ~issiso(sys)
        error('systemToRawQP requires a SISO system.');
    end

    sys = ss(sys);


    % ============================================================
    % Ensure OutputDelay does not enter internal reconstruction
    % ============================================================

    try
        sys.OutputDelay = 0;
    catch
    end


    % ============================================================
    % Internal-delay realization
    % ============================================================

    [A,B1,B2,C1,C2,D11,D12,D21,D22,E,tau] = ...
        getDelayModel(sys);

    tau = double(tau(:));


    % ============================================================
    % No internal delays
    % ============================================================

    if isempty(tau)

        G = tf(sys);

        [num,den] = tfdata(G,'v');

        num = cleanPolynomial(num,tol);
        den = cleanPolynomial(den,tol);

        R.numPoly  = {num};
        R.numDelay = 0;

        R.denPoly  = {den};
        R.denDelay = 0;

        R.internalDelays = [];

        return
    end


    % ============================================================
    % Symbolic variables
    %
    % one independent z_i for each internal delay:
    %
    %       z_i = exp(-tau_i*s)
    % ============================================================

    syms s

    nz = numel(tau);

    z = sym('z',[nz 1]);


    A   = sym(A);
    B1  = sym(B1);
    B2  = sym(B2);

    C1  = sym(C1);
    C2  = sym(C2);

    D11 = sym(D11);
    D12 = sym(D12);
    D21 = sym(D21);
    D22 = sym(D22);

    if isempty(E)
        E = eye(size(A));
    end

    E = sym(E);


    % ============================================================
    % Delay-free transfer blocks
    % ============================================================

    if isempty(A)

        H11 = D11;
        H12 = D12;
        H21 = D21;
        H22 = D22;

    else

        Rstate = inv(s*E-A);

        H11 = C1*Rstate*B1 + D11;
        H12 = C1*Rstate*B2 + D12;

        H21 = C2*Rstate*B1 + D21;
        H22 = C2*Rstate*B2 + D22;

    end


    % ============================================================
    % Internal delay operator
    % ============================================================

    Lambda = diag(z);

    Iz = sym(eye(nz));


    % ============================================================
    % Close internal-delay loop
    %
    % w = Lambda*zsignal
    %
    % zsignal = H21*u + H22*w
    %
    % therefore
    %
    % G = H11 + H12*(I-Lambda*H22)^(-1)*Lambda*H21
    % ============================================================

    Gsym = H11 + ...
        H12*((Iz-Lambda*H22)\(Lambda*H21));

    if numel(Gsym) ~= 1
        error('Internal realization did not produce a SISO expression.');
    end


    % ============================================================
    % Simplify symbolic rational expression
    % ============================================================

    Gsym = safeSimplifyFraction(Gsym);

    [Ns,Ds] = numden(Gsym);

    Ns = expand(Ns);
    Ds = expand(Ds);


    % ============================================================
    % Convert symbolic numerator / denominator into raw QP terms
    % ============================================================

    [numPoly,numDelay] = symbolicQPToTerms( ...
        Ns,z,tau,s,tol);

    [denPoly,denDelay] = symbolicQPToTerms( ...
        Ds,z,tau,s,tol);


    % ============================================================
    % Combine equal delays
    % ============================================================

    [numPoly,numDelay] = combineEqualDelays( ...
        numPoly,numDelay,tol);

    [denPoly,denDelay] = combineEqualDelays( ...
        denPoly,denDelay,tol);


    % ============================================================
    % Remove zero terms
    % ============================================================

    [numPoly,numDelay] = removeZeroTerms( ...
        numPoly,numDelay,tol);

    [denPoly,denDelay] = removeZeroTerms( ...
        denPoly,denDelay,tol);


    % ============================================================
    % Output
    % ============================================================

    R.numPoly  = numPoly;
    R.numDelay = numDelay;

    R.denPoly  = denPoly;
    R.denDelay = denDelay;

    R.internalDelays = tau;

end


% =================================================================
% SYMBOLIC QP -> RAW POLYNOMIAL/DELAY TERMS
% =================================================================
function [polyCell,delays] = symbolicQPToTerms( ...
    expr,z,tau,s,tol)

    expr = expand(expr);

    zrow = z(:).';

    [C,T] = coeffs(expr,zrow);

    C = C(:);
    T = T(:);

    polyCell = {};
    delays   = [];


    for j = 1:numel(C)

        cj = expand(C(j));
        tj = T(j);

        totalDelay = 0;


        % --------------------------------------------------------
        % Determine actual delay of delay monomial
        % --------------------------------------------------------

        for ii = 1:numel(z)

            pwr = symbolicDegree(tj,z(ii));

            totalDelay = totalDelay + ...
                double(pwr)*tau(ii);

        end


        % --------------------------------------------------------
        % Coefficient should depend only on s
        % --------------------------------------------------------

        otherVars = symvar(cj);

        for ii = 1:numel(otherVars)

            if ~isequal(otherVars(ii),s)

                error(['Coefficient still contains symbolic ', ...
                       'variable %s.'],char(otherVars(ii)));

            end

        end


        if isequal(cj,sym(0))
            continue
        end


        % --------------------------------------------------------
        % Convert polynomial in s to numeric coefficient vector
        % --------------------------------------------------------

        try

            p = double(sym2poly(cj));

        catch ME

            error(['Unable to convert symbolic coefficient to ', ...
                   'a polynomial in s:\n%s\n\n%s'], ...
                   char(cj),ME.message);

        end


        p = cleanPolynomial(p,tol);

        if all(abs(p)<tol)
            continue
        end


        polyCell{end+1,1} = p; %#ok<AGROW>
        delays(end+1,1) = totalDelay; %#ok<AGROW>

    end


    if isempty(polyCell)

        polyCell = {0};
        delays = 0;

    end

end


% =================================================================
% SYMBOLIC DEGREE
% =================================================================
function d = symbolicDegree(expr,x)

    try

        d = polynomialDegree(expr,x);

    catch

        d = feval(symengine,'degree',expr,x);

    end

    d = double(d);

    if isempty(d) || ~isfinite(d)
        d = 0;
    end

    d = round(d);

end


% =================================================================
% SAFE FRACTION SIMPLIFICATION
% =================================================================
function f = safeSimplifyFraction(f)

    try

        f = simplifyFraction(f,'Expand',true);

    catch

        try

            f = simplifyFraction(f);

        catch

            f = simplify(f);

            [n,d] = numden(f);

            f = expand(n)/expand(d);

        end

    end

end


% =================================================================
% MULTIPLY TWO QUASIPOLYNOMIALS
% =================================================================
function [Cpoly,Cdelay] = multiplyQP( ...
    Apoly,Adelay,Bpoly,Bdelay,tol)

    Cpoly  = {};
    Cdelay = [];


    for i = 1:numel(Apoly)

        for j = 1:numel(Bpoly)

            p = conv(Apoly{i},Bpoly{j});

            tau = Adelay(i) + Bdelay(j);

            Cpoly{end+1,1} = p; %#ok<AGROW>
            Cdelay(end+1,1) = tau; %#ok<AGROW>

        end

    end


    [Cpoly,Cdelay] = combineEqualDelays( ...
        Cpoly,Cdelay,tol);

end


% =================================================================
% COMBINE EQUAL-DELAY TERMS
% =================================================================
function [polyOut,delayOut] = combineEqualDelays( ...
    polyIn,delayIn,tol)

    if isempty(polyIn)

        polyOut = {};
        delayOut = [];

        return

    end


    delayIn = delayIn(:);

    [delayIn,idx] = sort(delayIn);

    polyIn = polyIn(idx);


    polyOut  = {};
    delayOut = [];


    for i = 1:numel(polyIn)

        tau = delayIn(i);

        p = cleanPolynomial(polyIn{i},tol);


        if isempty(delayOut)

            delayOut(1,1) = tau;
            polyOut{1,1} = p;

            continue

        end


        if abs(tau-delayOut(end)) <= ...
                tol*max([1,abs(tau),abs(delayOut(end))])

            polyOut{end} = addPolynomials( ...
                polyOut{end},p);

        else

            delayOut(end+1,1) = tau; %#ok<AGROW>
            polyOut{end+1,1} = p; %#ok<AGROW>

        end

    end


    [polyOut,delayOut] = removeZeroTerms( ...
        polyOut,delayOut,tol);

end


% =================================================================
% REMOVE ZERO TERMS
% =================================================================
function [polyOut,delayOut] = removeZeroTerms( ...
    polyIn,delayIn,tol)

    polyOut  = {};
    delayOut = [];


    for i = 1:numel(polyIn)

        p = cleanPolynomial(polyIn{i},tol);

        if all(abs(p)<=tol)
            continue
        end

        polyOut{end+1,1} = p; %#ok<AGROW>
        delayOut(end+1,1) = delayIn(i); %#ok<AGROW>

    end

end


% =================================================================
% ADD POLYNOMIALS
% =================================================================
function p = addPolynomials(a,b)

    a = a(:).';
    b = b(:).';

    L = max(numel(a),numel(b));

    a = [zeros(1,L-numel(a)) a];
    b = [zeros(1,L-numel(b)) b];

    p = a+b;

end


% =================================================================
% NUMERICAL GCD OF INTERNAL DELAYS
% =================================================================
function h = delayGCD(tau,tol)

    tau = tau(:);

    tau(abs(tau)<tol) = 0;

    tp = tau(tau>tol);


    % Pure polynomial
    if isempty(tp)

        h = 0;
        return

    end


    digits = max(0,ceil(-log10(tol)));
    digits = min(digits,9);

    scale = 10^digits;

    ti = round(tp*scale);

    g = abs(ti(1));


    for i = 2:numel(ti)

        g = gcd(g,abs(ti(i)));

    end


    h = double(g)/scale;


    if h <= tol

        error('Unable to determine a useful internal delay GCD.');

    end


    % Verify
    k = round(tp/h);

    err = abs(tp-k*h);

    if any(err > 100*tol.*max(1,abs(tp)))

        error(['Internal delays are not commensurate. ', ...
               'Maximum residual = %.4e.'],max(err));

    end

end


% =================================================================
% PUT QP ON REGULAR INTERNAL-DELAY GRID
% =================================================================
function [M,tauGrid,k] = qpToGrid( ...
    polyCell,delays,h,tol)

    delays = delays(:);


    % ============================================================
    % Pure polynomial
    % ============================================================

    if h == 0

        tauGrid = 0;

        k = zeros(size(delays));

        p = 0;

        for i = 1:numel(polyCell)

            p = addPolynomials(p,polyCell{i});

        end

        M = cleanPolynomial(p,tol);

        return

    end


    % ============================================================
    % Integer delay multipliers
    % ============================================================

    k = round(delays/h);

    err = abs(delays-k*h);


    if any(err > 100*tol.*max(1,abs(delays)))

        error(['Delay is not an integer multiple of ', ...
               'h = %.15g.'],h);

    end


    Kmax = max(k);

    tauGrid = (0:Kmax).' * h;


    % ============================================================
    % Maximum polynomial degree
    % ============================================================

    Lmax = 1;

    for i = 1:numel(polyCell)

        Lmax = max(Lmax,numel(polyCell{i}));

    end


    M = zeros(Kmax+1,Lmax);


    % ============================================================
    % Fill rows
    % ============================================================

    for i = 1:numel(polyCell)

        p = polyCell{i}(:).';

        p = [zeros(1,Lmax-numel(p)) p];

        row = k(i)+1;

        M(row,:) = M(row,:) + p;

    end


    M(abs(M)<tol) = 0;

end


% =================================================================
% CLEAN POLYNOMIAL
% =================================================================
function p = cleanPolynomial(p,tol)

    p = double(p(:).');

    if isempty(p)

        p = 0;
        return

    end


    if max(abs(imag(p))) <= ...
            100*tol*max(1,max(abs(p)))

        p = real(p);

    end


    p(abs(p)<tol) = 0;


    idx = find(abs(p)>tol,1,'first');


    if isempty(idx)

        p = 0;

    else

        p = p(idx:end);

    end

end


% =================================================================
% FIRST NONZERO COEFFICIENT
% =================================================================
function c = firstNonzeroCoefficient(p,tol)

    idx = find(abs(p)>tol,1,'first');

    if isempty(idx)

        c = [];

    else

        c = p(idx);

    end

end