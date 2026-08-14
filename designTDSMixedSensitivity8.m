function [K,info] = designTDSMixedSensitivity8(C0,varargin)
% DESIGNTDSMIXEDSENSITIVITY8_GRID15
%
% Frequency-grid HANSO optimization for the exact delayed plant.
%
% Controller structure
% --------------------
%   K(s) = Q(s)/(s+1e-5)
%
% where Q is 7th-order and biproper:
%
%   Q(s) = (b7*s^7+...+b0)/(s^7+a6*s^6+...+a0).
%
% Hence K is order 8, relative degree 1, and contains a fixed pole -1e-5.
% There are 15 free parameters:
%   7 denominator coefficients + 8 numerator coefficients.
%
% Initialization
% --------------
% The closest REAL pole of C0 to -1e-5 is replaced by exactly -1e-5:
%
%   Kinit = C0*(s-pclose)/(s+1e-5),
%   Q0    = (s+1e-5)*Kinit.
%
% Objective
% ---------
% No tds_hinfnorm call is used during optimization.
%
% On a dense frequency grid:
%
%   S = 1/(1+P*K),    T = P*K/(1+P*K)
%
% and HANSO minimizes
%
%   max_w sqrt(|W1*S|^2 + |W2*T|^2),
%
% using the exact analytical delayed frequency response P(jw).
%
% Stability
% ---------
% TDS-CONTROL is used only to check characteristic roots of the exact
% delayed closed loop with tds_create_cl + tds_roots.  Unstable trial
% controllers receive a large penalty.  This avoids the repeated RCOND
% warnings produced by tds_hinfnorm.
%
% OPTIONS
% -------
%   'Wmin'            : 1e-3 rad/s
%   'Wmax'            : 1e4  rad/s
%   'Nw'              : 4000 logarithmic samples (+ w=0)
%   'StabilityMargin' : 0
%   'RootSearchLeft'  : -1
%   'Quiet'           : true
%
% OUTPUTS
% -------
%   K    : optimized order-8 controller
%   info : diagnostics, grid performance and stability information
%
% NOTE
% ----
% HANSO 3.0 requires pars.fgname to name an m-file. To keep this as a
% single user-facing m-file, a tiny callback m-file is generated
% automatically under tempdir while HANSO runs.

    ip = inputParser;
    ip.addParameter('Wmin',1e-3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    ip.addParameter('Wmax',1e4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    ip.addParameter('Nw',4000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
    ip.addParameter('StabilityMargin',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    ip.addParameter('RootSearchLeft',-1,@(x)isnumeric(x)&&isscalar(x)&&x<0);
    ip.addParameter('Quiet',true,@(x)islogical(x)&&isscalar(x));
    ip.parse(varargin{:});
    o = ip.Results;

    pfix = 1e-5;
    s = tf('s');

    %% Initial controller -> fixed-pole family
    C0 = minreal(tf(C0));

    if ~issiso(C0)
        error('C0 must be SISO.');
    end

    p0 = pole(C0);
    ir = find(abs(imag(p0)) < 1e-8);
    if isempty(ir)
        error('C0 has no real pole to replace by -1e-5.');
    end

    [~,jj] = min(abs(real(p0(ir)) + pfix));
    pclose = real(p0(ir(jj)));

    Kinit = minreal(C0*(s-pclose)/(s+pfix));
    Q0 = minreal((s+pfix)*Kinit);

    [nq0,dq0] = tfdata(Q0,'v');
    nq0 = nq0/dq0(1);
    dq0 = dq0/dq0(1);

    % Remove only numerical leading zeros.
    kn = find(abs(nq0)>1e-12*max(1,norm(nq0,inf)),1,'first');
    nq0 = nq0(kn:end);

    if numel(dq0) ~= 8 || numel(nq0) ~= 8
        error('Q0 must be 7th-order and biproper.');
    end

    % Physical parameter vector:
    %   theta = [a6 ... a0  b7 ... b0]
    theta0 = [dq0(2:end), nq0].';

    % Dimensionless HANSO coordinates to reduce coefficient-scale problems.
    denScale = max(abs(theta0(1:7)), max(abs(theta0(1:7)))*1e-4 + 1e-12);
    numScale = max(abs(theta0(8:15)), max(abs(theta0(8:15)))*1e-4 + 1e-12);
    scale = [denScale; numScale];
    z0 = theta0./scale;

    %% Exact delayed plant for TDS-CONTROL stability checks
    tau = [0 0.2 0.3 1.5 2];

    Pnum = [ ...
         0   0;
         0   5;
        -1   0;
         0   0;
         0   0];

    Pden = [ ...
         3   3.5   0.5;
         0   0     0;
         0   0     0;
         0   2     7;
         0   1    -1];

    Ptds = tds_create_tf(Pnum,Pden,tau);

    %% Exact frequency grid
    w = [0, logspace(log10(o.Wmin),log10(o.Wmax),round(o.Nw))];
    jw = 1i*w;

    % Exact delayed P(jw), no approximation of any delay.
    Pw = exp(-0.2*jw).*(5 - jw.*exp(-0.1*jw)) ./ ...
        ( (jw+1).*(3*jw+0.5) ...
        + (2*jw+7).*exp(-1.5*jw) ...
        + (jw-1).*exp(-2*jw) );

    W1w = 0.5./(jw+pfix);
    W2w = 0.01*jw.^2 + 0.3*jw + 0.2;

    %% Global context for HANSO callback
    global TDS_GRID15_CTX
    TDS_GRID15_CTX = struct;
    TDS_GRID15_CTX.w = w;
    TDS_GRID15_CTX.jw = jw;
    TDS_GRID15_CTX.Pw = Pw;
    TDS_GRID15_CTX.W1w = W1w;
    TDS_GRID15_CTX.W2w = W2w;
    TDS_GRID15_CTX.Ptds = Ptds;
    TDS_GRID15_CTX.pfix = pfix;
    TDS_GRID15_CTX.scale = scale;
    TDS_GRID15_CTX.margin = o.StabilityMargin;
    TDS_GRID15_CTX.rootLeft = o.RootSearchLeft;
    TDS_GRID15_CTX.quiet = o.Quiet;

    %% Verify structured initial controller
    alpha0 = localClosedLoopAlpha(Ptds,Kinit,o.RootSearchLeft,o.Quiet);
    [g0,k0] = localGridObjective(theta0,w,jw,Pw,W1w,W2w,pfix);

    fprintf('Closest real pole replaced: %.12e -> %.12e\n',pclose,-pfix);
    fprintf('Initial delayed-loop alpha = %.12e\n',alpha0);
    fprintf('Initial grid gamma         = %.12e at %.12e rad/s\n',g0,w(k0));

    if ~isfinite(alpha0) || alpha0 >= -o.StabilityMargin
        clear global TDS_GRID15_CTX
        error(['The pole-replaced initial controller is not stable with the ' ...
               'requested delayed closed-loop stability margin.']);
    end

    %% HANSO 3.0
    [fgname,cleanupObj] = localEnsureHansoCallback(); %#ok<NASGU>

    pars = struct;
    pars.fgname = fgname;
    pars.nvar = 15;

    options = struct;
    options.x0 = z0;
    options.prtlevel = 1;

    [zBest,fBest] = hanso(pars,options);

    thetaBest = scale.*zBest(:);
    K = localThetaToK(thetaBest,pfix);

    %% Final diagnostics
    alpha = localClosedLoopAlpha(Ptds,K,o.RootSearchLeft,o.Quiet);
    [gamma,kPeak,perf] = localGridObjective(thetaBest,w,jw,Pw,W1w,W2w,pfix);

    info = struct;
    info.gammaGrid = gamma;
    info.hansoObjective = fBest;
    info.peakFrequency = w(kPeak);
    info.frequency = w;
    info.performance = perf;
    info.closedLoopAlpha = alpha;
    info.initialGammaGrid = g0;
    info.initialPeakFrequency = w(k0);
    info.initialClosedLoopAlpha = alpha0;
    info.Kinit = Kinit;
    info.Q0 = Q0;
    info.replacedPole = pclose;
    info.fixedPole = -pfix;
    info.parameters = thetaBest;
    info.parameterScale = scale;
    info.nParameters = 15;
    info.Ptds = Ptds;
    info.W1 = 0.5/(s+pfix);
    info.W2 = 0.01*s^2+0.3*s+0.2;

    clear global TDS_GRID15_CTX
end


function [gamma,kPeak,perf] = localGridObjective(theta,w,jw,Pw,W1w,W2w,pfix) %#ok<INUSD>
    a = theta(1:7).';
    b = theta(8:15).';

    D = polyval([1 a],jw);
    N = polyval(b,jw);

    Qw = N./D;
    Kw = Qw./(jw+pfix);

    L = Pw.*Kw;
    S = 1./(1+L);
    T = L.*S;

    perf = sqrt(abs(W1w.*S).^2 + abs(W2w.*T).^2);
    [gamma,kPeak] = max(perf);
end


function K = localThetaToK(theta,pfix)
    a = theta(1:7).';
    b = theta(8:15).';

    denQ = [1 a];
    denK = conv([1 pfix],denQ);

    K = ss(tf(b,denK));
end


function alpha = localClosedLoopAlpha(P,K,left,quiet)
% Build the dynamic P-K closed loop explicitly.
%
% tds_create_cl in this TDS-CONTROL release treats a MATLAB ss object as a
% static D11 object on this path.  Since K is strictly proper, construct the
% augmented delayed dynamics directly:
%
%   xPdot = A_P xP + B_P C_K xK
%   xKdot = A_K xK - B_K(C_P xP + D_P C_K xK)

    CL = localBuildDynamicCL(P,K);
    ro = tds_roots_options('quiet',logical(quiet));

    bounds = unique([left, -5, -20, -100],'stable');
    alpha = NaN;

    for r = bounds
        if r >= 0
            continue;
        end
        try
            rr = tds_roots(CL,r,ro);
            lam = localExtractRoots(rr);
            if ~isempty(lam)
                alpha = max(real(lam));
                return;
            end
        catch
        end
    end
end


function CL = localBuildDynamicCL(P,K)
    K = ss(K);
    [Ak,Bk,Ck,Dk] = ssdata(K);

    if norm(Dk,'fro') > 1e-10
        error('The stability helper currently requires strictly proper K (D=0).');
    end

    Ap = localAsCell(P.A);   hAp = localRow(P.hA);
    Bp = localAsCell(P.B1);  hBp = localRow(P.hB1);
    Cp = localAsCell(P.C1);  hCp = localRow(P.hC1);

    if isprop(P,'D11') && ~isempty(P.D11)
        Dp = localAsCell(P.D11);
        hDp = localRow(P.hD11);
    else
        Dp = {};
        hDp = [];
    end

    np = size(Ap{1},1);
    nk = size(Ak,1);
    n = np+nk;

    ip = 1:np;
    ik = np+(1:nk);

    Acl = {};
    hAcl = [];

    % Plant autonomous terms.
    for j=1:numel(Ap)
        M=zeros(n);
        M(ip,ip)=Ap{j};
        [Acl,hAcl]=localAddTerm(Acl,hAcl,M,hAp(j));
    end

    % Controller autonomous term.
    M=zeros(n);
    M(ik,ik)=Ak;
    [Acl,hAcl]=localAddTerm(Acl,hAcl,M,0);

    % Plant input u = Ck*xk.
    for j=1:numel(Bp)
        M=zeros(n);
        M(ip,ik)=Bp{j}*Ck;
        [Acl,hAcl]=localAddTerm(Acl,hAcl,M,hBp(j));
    end

    % Controller input y = -P*u.
    for j=1:numel(Cp)
        M=zeros(n);
        M(ik,ip)=-Bk*Cp{j};
        [Acl,hAcl]=localAddTerm(Acl,hAcl,M,hCp(j));
    end

    % Plant direct delayed terms: -Bk*Dp*Ck*xk.
    for j=1:numel(Dp)
        M=zeros(n);
        M(ik,ik)=-Bk*Dp{j}*Ck;
        [Acl,hAcl]=localAddTerm(Acl,hAcl,M,hDp(j));
    end

    CL = tds_create(Acl,hAcl);
end


function c = localAsCell(x)
    if isempty(x)
        c={};
    elseif iscell(x)
        c=x;
    else
        c={x};
    end
end


function x = localRow(x)
    x=x(:).';
end


function [C,h]=localAddTerm(C,h,M,t)
    if isempty(M) || norm(M,'fro')<1e-14
        return
    end

    if isempty(h)
        C={M};
        h=t;
        return
    end

    j=find(abs(h-t)<=1e-12*max(1,abs(t)),1);
    if isempty(j)
        C{end+1}=M;
        h(end+1)=t;
    else
        C{j}=C{j}+M;
    end
end


function lam = localExtractRoots(rr)
    if isnumeric(rr)
        if isreal(rr) && size(rr,2)==2 && size(rr,1)>1
            lam = rr(:,1)+1i*rr(:,2);
        else
            lam = rr(:);
        end
        return;
    end

    names = {'roots','root','lambda','values','eig','eigenvalues'};
    for k=1:numel(names)
        f=names{k};
        try
            if isstruct(rr) && isfield(rr,f)
                lam=rr.(f); lam=lam(:); return;
            elseif isobject(rr) && isprop(rr,f)
                lam=rr.(f); lam=lam(:); return;
            end
        catch
        end
    end
    lam=[];
end


function [fgname,cleanupObj] = localEnsureHansoCallback()
% HANSO 3.0 requires pars.fgname to name a visible m-file. Generate the
% callback under tempdir so the user only needs this one source file.

    fgname = 'tdsGrid15FG_runtime';
    d = fullfile(tempdir,'tds_grid15_hanso');
    if ~exist(d,'dir')
        mkdir(d);
    end

    f = fullfile(d,[fgname '.m']);
    src = localCallbackSource();

    rewrite = true;
    if exist(f,'file')
        try
            rewrite = ~strcmp(fileread(f),src);
        catch
        end
    end

    if rewrite
        fid = fopen(f,'w');
        if fid<0
            error('Could not create temporary HANSO callback.');
        end
        fwrite(fid,src,'char');
        fclose(fid);
    end

    addpath(d,'-begin');
    cleanupObj = onCleanup(@() localRmPath(d));
end


function localRmPath(d)
    try
        rmpath(d);
    catch
    end
end


function src = localCallbackSource()
% Runtime HANSO callback. The stable-region objective/gradient is computed
% analytically on the grid. tds_roots is used only for the stability test.

L = {
'function [f,g] = tdsGrid15FG_runtime(z,varargin)'
'global TDS_GRID15_CTX'
'c=TDS_GRID15_CTX;'
'z=z(:);'
'theta=c.scale.*z;'
''
'% Exact delayed closed-loop stability check'
'K=thetaToK(theta,c.pfix);'
'alpha=closedLoopAlpha(c.Ptds,K,c.rootLeft,c.quiet);'
'if ~isfinite(alpha)'
'    f=1e8; g=zeros(15,1); return'
'end'
'if alpha >= -c.margin'
'    v=alpha+c.margin;'
'    f=1e6+1e5*v+1e4*v^2;'
'    g=zeros(15,1);'
'    return'
'end'
''
'% Grid objective'
'a=theta(1:7).'';'
'b=theta(8:15).'';'
'D=polyval([1 a],c.jw);'
'N=polyval(b,c.jw);'
'Q=N./D;'
'Kf=Q./(c.jw+c.pfix);'
'Lw=c.Pw.*Kf;'
'S=1./(1+Lw);'
'T=Lw.*S;'
'perf=sqrt(abs(c.W1w.*S).^2+abs(c.W2w.*T).^2);'
'[f,k]=max(perf);'
''
'if nargout<2, return; end'
''
'% Analytic active-peak gradient'
'sk=c.jw(k);'
'Pk=c.Pw(k);'
'S0=S(k); T0=T(k);'
'D0=D(k); N0=N(k);'
'phi=max(f,1e-15);'
''
'% d phi / d L for real parameters'
'coef=(S0^2/phi)*(-abs(c.W1w(k))^2*conj(S0)+abs(c.W2w(k))^2*conj(T0));'
''
'bd=sk.^(6:-1:0);'
'bn=sk.^(7:-1:0);'
'dQda=-(N0/(D0^2))*bd;'
'dQdb=(1/D0)*bn;'
'dK=[dQda dQdb]/(sk+c.pfix);'
'dL=Pk*dK;'
'gtheta=real(coef*dL).'';'
'g=c.scale.*gtheta;'
'end'
''
'function K=thetaToK(theta,pfix)'
'a=theta(1:7).''; b=theta(8:15).'';'
'K=ss(tf(b,conv([1 pfix],[1 a])));'
'end'
''
'function alpha=closedLoopAlpha(P,K,left,quiet)'
'CL=buildDynamicCL(P,K);'
'ro=tds_roots_options(''quiet'',logical(quiet));'
'bounds=unique([left -5 -20 -100],''stable'');'
'alpha=NaN;'
'for r=bounds'
'    if r>=0, continue; end'
'    try'
'        rr=tds_roots(CL,r,ro);'
'        lam=extractRoots(rr);'
'        if ~isempty(lam), alpha=max(real(lam)); return; end'
'    catch'
'    end'
'end'
'end'
''
'function CL=buildDynamicCL(P,K)'
'K=ss(K); [Ak,Bk,Ck,Dk]=ssdata(K);'
'if norm(Dk,''fro'')>1e-10, error(''K must be strictly proper.''); end'
'Ap=cc(P.A); hAp=rv(P.hA);'
'Bp=cc(P.B1); hBp=rv(P.hB1);'
'Cp=cc(P.C1); hCp=rv(P.hC1);'
'if isprop(P,''D11'') && ~isempty(P.D11)'
'    Dp=cc(P.D11); hDp=rv(P.hD11);'
'else'
'    Dp={}; hDp=[];'
'end'
'np=size(Ap{1},1); nk=size(Ak,1); n=np+nk;'
'ip=1:np; ik=np+(1:nk);'
'Acl={}; hAcl=[];'
'for j=1:numel(Ap)'
'    M=zeros(n); M(ip,ip)=Ap{j}; [Acl,hAcl]=at(Acl,hAcl,M,hAp(j));'
'end'
'M=zeros(n); M(ik,ik)=Ak; [Acl,hAcl]=at(Acl,hAcl,M,0);'
'for j=1:numel(Bp)'
'    M=zeros(n); M(ip,ik)=Bp{j}*Ck; [Acl,hAcl]=at(Acl,hAcl,M,hBp(j));'
'end'
'for j=1:numel(Cp)'
'    M=zeros(n); M(ik,ip)=-Bk*Cp{j}; [Acl,hAcl]=at(Acl,hAcl,M,hCp(j));'
'end'
'for j=1:numel(Dp)'
'    M=zeros(n); M(ik,ik)=-Bk*Dp{j}*Ck; [Acl,hAcl]=at(Acl,hAcl,M,hDp(j));'
'end'
'CL=tds_create(Acl,hAcl);'
'end'
''
'function c=cc(x)'
'if isempty(x), c={}; elseif iscell(x), c=x; else, c={x}; end'
'end'
'function x=rv(x), x=x(:).''; end'
'function [C,h]=at(C,h,M,t)'
'if isempty(M)||norm(M,''fro'')<1e-14, return; end'
'if isempty(h), C={M}; h=t; return; end'
'j=find(abs(h-t)<=1e-12*max(1,abs(t)),1);'
'if isempty(j), C{end+1}=M; h(end+1)=t; else, C{j}=C{j}+M; end'
'end'
''
'function lam=extractRoots(rr)'
'if isnumeric(rr)'
'    if isreal(rr) && size(rr,2)==2 && size(rr,1)>1'
'        lam=rr(:,1)+1i*rr(:,2);'
'    else'
'        lam=rr(:);'
'    end'
'    return'
'end'
'names={''roots'',''root'',''lambda'',''values'',''eig'',''eigenvalues''};'
'for k=1:numel(names)'
'    q=names{k};'
'    try'
'        if isstruct(rr) && isfield(rr,q)'
'            lam=rr.(q); lam=lam(:); return'
'        elseif isobject(rr) && isprop(rr,q)'
'            lam=rr.(q); lam=lam(:); return'
'        end'
'    catch'
'    end'
'end'
'lam=[];'
'end'
};
src=sprintf('%s\n',L{:});
end
