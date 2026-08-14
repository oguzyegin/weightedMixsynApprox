function [Cm,gammaBest,info] = tfirkaMixedSensitivity(Copt,P,Pa,W1,W2,m,opts)
%TFIRKAMIXEDSENSITIVITY
% Optimize TF-IRKA interpolation points for an exactly order-m controller
% that stabilizes Pa and minimizes the true mixed-sensitivity peak.
%
%   [Cm,gammaBest,info] = tfirkaMixedSensitivity( ...
%       Copt,P,Pa,W1,W2,m)
%
% INPUTS
%   Copt : symbolic scalar expression in Laplace variable s. This is the
%          high/infinite-dimensional controller to be approximated.
%
%   P    : plant used for the mixed-sensitivity objective. May be tf/zpk/ss/
%          dss or FRD. If FRD is supplied, the objective is evaluated on its
%          frequency grid.
%
%   Pa   : finite-dimensional SISO plant used for the hard closed-loop
%          stability test.
%
%   W1   : sensitivity weight. tf/zpk/ss/dss or FRD.
%   W2   : complementary-sensitivity weight. tf/zpk/ss/dss or FRD.
%
%   m    : EXACT desired controller order.
%
% OPTIONS
%   .Frequency
%       Frequency grid for the mixed-sensitivity objective. If empty:
%         - use P.Frequency when P is FRD;
%         - otherwise use logspace(-5,4,2500).
%
%   .InitialPoints
%       Optional complete m-point complex interpolation set. Must be
%       conjugation-closed for a real controller.
%
%   .UseUnstablePoleInitialization   default true
%   .UseSensitivityPeakInitialization default true
%   .InitialRealPart                default 1
%       Default positive real part of automatically initialized TF-IRKA
%       interpolation points.
%
%   .StabilityMargin                default 1e-7
%       Require max(real(pole(feedback(Pa*Cm,1)))) < -StabilityMargin.
%
%   .MaxIterations                  default 150
%   .MaxFunctionEvaluations         default 5000
%   .Display                        'iter'|'off'|'final', default 'iter'
%   .Optimizer                      'fminsearch'|'particleswarm', default
%                                   'fminsearch'. particleswarm requires
%                                   Global Optimization Toolbox.
%
%   .PointRealLowerBound            default 1e-5
%   .PointRealUpperBound            default 1e4
%   .PointImagUpperBound            default 1e5
%   .FailurePenalty                 default 1e8
%   .UnstablePenaltySlope           default 1e6
%   .Regularization                 default 1e-6
%
% OUTPUTS
%   Cm        : best exactly order-m TF-IRKA controller found.
%   gammaBest : best mixed-sensitivity peak.
%   info      : optimization history, interpolation points, stability data.
%
% MIXED-SENSITIVITY OBJECTIVE
%   S = 1/(1+P*Cm)
%   T = P*Cm/(1+P*Cm)
%
%   gamma = max_w sqrt(|W1*S|^2 + |W2*T|^2)
%
% CONTROLLER CONSTRUCTION
%   Every candidate controller is constructed ONLY through the Hermite/
%   Loewner TF-IRKA realization from Copt(s),Copt'(s) at the current
%   interpolation points.
%
% REAL CONTROLLER PARAMETERIZATION
%   For even m=2r:
%       sigma_k = exp(alpha_k) + j*exp(beta_k),  k=1,...,r
%       plus conjugates.
%
%   For odd m=2r+1:
%       one positive real point exp(alpha_0), plus r conjugate pairs.
%
%   Thus optimization variables are unconstrained real logarithmic
%   coordinates and the interpolation points remain in the RHP.

arguments
    Copt sym
    P
    Pa
    W1
    W2
    m (1,1) double {mustBeInteger,mustBePositive}
    opts.Frequency double = []
    opts.InitialPoints double = []
    opts.UseUnstablePoleInitialization (1,1) logical = true
    opts.UseSensitivityPeakInitialization (1,1) logical = true
    opts.InitialRealPart (1,1) double {mustBePositive} = 1
    opts.StabilityMargin (1,1) double {mustBeNonnegative} = 1e-6
    opts.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 2000
    opts.MaxFunctionEvaluations (1,1) double {mustBeInteger,mustBePositive} = 5000
    opts.Display char {mustBeMember(opts.Display,{'iter','off','final'})} = 'iter'
    opts.Optimizer char {mustBeMember(opts.Optimizer,{'fminsearch','particleswarm'})} = 'fminsearch'
    opts.PointRealLowerBound (1,1) double {mustBePositive} = 1e-5
    opts.PointRealUpperBound (1,1) double {mustBePositive} = 1e4
    opts.PointImagUpperBound (1,1) double {mustBePositive} = 1e5
    opts.FailurePenalty (1,1) double {mustBePositive} = 1e8
    opts.UnstablePenaltySlope (1,1) double {mustBePositive} = 1e6
    opts.Regularization (1,1) double {mustBeNonnegative} = 1e-6
end
if isempty(opts.InitialPoints)
    pu = pole(Pa);
    pu = pu(real(pu) > 0);

    opts.InitialPoints = [ ...
        pu.', ...
        0.02+0.7j, 0.02-0.7j, ...
        0.1-4.2j, 0.1+4.2j, ...
        1-15j,   1+15j ];
end
%% Validate Pa
if ~(isa(Pa,'tf') || isa(Pa,'zpk') || isa(Pa,'ss') || isa(Pa,'dss'))
    error('Pa must be a finite-dimensional tf/zpk/ss/dss plant.');
end
if size(Pa,1)~=1 || size(Pa,2)~=1 || Pa.Ts~=0
    error('Pa must be continuous-time SISO.');
end

%% Symbolic evaluators
sv = symvar(Copt);
if isempty(sv)
    error('Copt must depend on a symbolic Laplace variable.');
end
s = sv(1);

CoptPrime = diff(Copt,s);
Cfun = matlabFunction(Copt,'Vars',s);
CpFun = matlabFunction(CoptPrime,'Vars',s);

%% Frequency grid
if isempty(opts.Frequency)
    if isa(P,'frd')
        w = P.Frequency(:);
    else
        w = logspace(-5,4,2500).';
    end
else
    w = opts.Frequency(:);
end

w = unique(w(isfinite(w) & w>0));
if isempty(w)
    error('Frequency grid must contain positive finite frequencies.');
end

%% Evaluate fixed plant/weights on objective grid
Pjw  = responseOnGrid(P,w);
W1jw = responseOnGrid(W1,w);
W2jw = responseOnGrid(W2,w);

%% Optional Copt FRD on grid for initialization only
Coptjw = evalSymbolic(Cfun,1i*w);
Loptjw = Pjw .* Coptjw;
Soptjw = 1 ./ (1 + Loptjw);

%% Initial interpolation set
if isempty(opts.InitialPoints)
    spt0 = automaticInitialPoints( ...
        Pa,w,Soptjw,m,opts.InitialRealPart, ...
        opts.UseUnstablePoleInitialization, ...
        opts.UseSensitivityPeakInitialization);
else
    spt0 = opts.InitialPoints(:);
    if numel(spt0) ~= m
        error('InitialPoints must contain exactly m points.');
    end
    spt0 = sanitizeInitialPointSet(spt0,m,opts);
end

x0 = pointsToParameters(spt0,m,opts);

historyX = {};
historyGamma = [];
historyStable = [];
historyMaxRealPole = [];
historyPoints = {};

best.gamma = Inf;
best.Cm = [];
best.points = [];
best.maxRealPole = Inf;
best.closedLoopPoles = [];

evalCount = 0;

%% Objective wrapper
    function J = outerCost(x)
        evalCount = evalCount + 1;

        try
            spt = parametersToPoints(x,m,opts);

            [Cmc,buildInfo] = tfirkaAtPoints(Cfun,CpFun,spt,m);

            % Exact-order guard
            if buildInfo.effectiveOrder ~= m
                J = opts.FailurePenalty + 1e5*(m-buildInfo.effectiveOrder)^2;
                record(J,false,Inf,spt);
                return
            end

            % Closed-loop stability on Pa
            CL = feedback(Pa*Cmc,1);
            pcl = pole(CL);

            if isempty(pcl)
                maxRealPole = -Inf;
                stable = true;
            else
                maxRealPole = max(real(pcl));
                stable = all(isfinite(pcl)) && ...
                         maxRealPole < -opts.StabilityMargin;
            end

            if ~stable
                violation = max(0,maxRealPole + opts.StabilityMargin);
                J = opts.FailurePenalty + ...
                    opts.UnstablePenaltySlope*violation^2;
                record(J,false,maxRealPole,spt);
                return
            end

            % True mixed-sensitivity objective
            Cjw = squeeze(freqresp(Cmc,w)); Cjw=Cjw(:);
            L = Pjw .* Cjw;
            den = 1 + L;

            if any(abs(den)<1e-12) || any(~isfinite(den))
                J = opts.FailurePenalty;
                record(J,false,maxRealPole,spt);
                return
            end

            Sresp = 1 ./ den;
            Tresp = L ./ den;

            g = sqrt(abs(W1jw.*Sresp).^2 + abs(W2jw.*Tresp).^2);
            gamma = max(g);

            % Gentle regularization discourages runaway point movement.
            reg = opts.Regularization*norm(x-x0)^2;
            J = gamma + reg;

            record(gamma,true,maxRealPole,spt);

            if gamma < best.gamma
                best.gamma = gamma;
                best.Cm = Cmc;
                best.points = spt;
                best.maxRealPole = maxRealPole;
                best.closedLoopPoles = pcl;
                best.buildInfo = buildInfo;
                best.gCurve = g;
                best.S = Sresp;
                best.T = Tresp;
            end

        catch
            J = opts.FailurePenalty;
            try
                spt = parametersToPoints(x,m,opts);
            catch
                spt = nan(m,1);
            end
            record(J,false,Inf,spt);
        end
    end

    function record(gam,stab,maxrp,spt)
        historyGamma(end+1,1) = gam;
        historyStable(end+1,1) = stab;
        historyMaxRealPole(end+1,1) = maxrp;
        historyPoints{end+1,1} = spt;
    end

%% Optimize interpolation parameters
switch opts.Optimizer
    case 'fminsearch'
        fopts = optimset( ...
            'Display',opts.Display, ...
            'MaxIter',opts.MaxIterations, ...
            'MaxFunEvals',opts.MaxFunctionEvaluations, ...
            'TolX',1e-5, ...
            'TolFun',1e-5);

        [xopt,fopt,exitflag,output] = fminsearch(@outerCost,x0,fopts);

    case 'particleswarm'
        nvar = numel(x0);

        % bounds are in log coordinates
        [lb,ub] = parameterBounds(m,opts);

        psopts = optimoptions('particleswarm', ...
            'Display',opts.Display, ...
            'MaxIterations',opts.MaxIterations, ...
            'MaxStallIterations',30, ...
            'FunctionTolerance',1e-5);

        [xopt,fopt,exitflag,output] = particleswarm( ...
            @outerCost,nvar,lb,ub,psopts);
end

%% Evaluate final optimizer point too
outerCost(xopt);

if isempty(best.Cm)
    error('tfirkaMixedSensitivity:NoStableController', ...
        ['No stabilizing order-%d TF-IRKA controller was found. ', ...
         'Try different initial points or a global optimizer.'],m);
end

Cm = best.Cm;
gammaBest = best.gamma;

%% Final information
info = struct;
info.gamma = gammaBest;
info.controllerOrder = order(ss(Cm));
info.interpolationPoints = best.points;
info.initialInterpolationPoints = spt0;
info.optimizedParameters = xopt;
info.optimizerObjective = fopt;
info.exitflag = exitflag;
info.output = output;
info.evaluationCount = evalCount;

info.closedLoopStable = true;
info.closedLoopPoles = best.closedLoopPoles;
info.maxRealClosedLoopPole = best.maxRealPole;

info.frequency = w;
info.gammaCurve = best.gCurve;
info.S = best.S;
info.T = best.T;
info.Pjw = Pjw;
info.W1jw = W1jw;
info.W2jw = W2jw;

info.historyGamma = historyGamma;
info.historyStable = logical(historyStable);
info.historyMaxRealPole = historyMaxRealPole;
info.historyPoints = historyPoints;

info.tfIRKABuildInfo = best.buildInfo;
end


function [Cm,info] = tfirkaAtPoints(Cfun,CpFun,spt,m)
% Construct the Hermite/Loewner TF-IRKA interpolant at a fixed point set.

spt = spt(:);

yi = evalSymbolic(Cfun,spt);
yp = evalSymbolic(CpFun,spt);

if any(~isfinite(yi)) || any(~isfinite(yp))
    error('Nonfinite Copt or derivative evaluation.');
end

E = zeros(m);
A = zeros(m);

for i=1:m
    for j=1:m
        if i==j
            E(i,i) = -yp(i);
            A(i,i) = -(yi(i)+spt(i)*yp(i));
        else
            ds = spt(i)-spt(j);
            if abs(ds) < 1e-9*max([1,abs(spt(i)),abs(spt(j))])
                error('Interpolation points too close.');
            end
            E(i,j) = -(yi(i)-yi(j))/ds;
            A(i,j) = -(spt(i)*yi(i)-spt(j)*yi(j))/ds;
        end
    end
end

B = yi(:);
C = yi(:).';

% Realify via conjugate-pair basis
[Er,Ar,Br,Cr,resid] = realifyGeneral(E,A,B,C,spt);

if rcond(Er) < 1e-12
    Gd = dss(Ar,Br,Cr,0,Er);
    G = tf(Gd);
else
    As = Er\Ar;
    Bs = Er\Br;
    G = tf(ss(As,Bs,Cr,0));
end

[num,den] = tfdata(G,'v');
num = real(num);
den = real(den);

% normalize
if den(1)~=0
    num = num/den(1);
    den = den/den(1);
end

Cm = minreal(tf(num,den),1e-9);

info = struct;
info.E = Er; info.A = Ar; info.B = Br; info.C = Cr;
info.realificationResidual = resid;
info.effectiveOrder = order(ss(Cm));
info.numerator = num;
info.denominator = den;
end


function spt = automaticInitialPoints(Pa,w,Sopt,m,realPart,usePu,useS)
% Build a conjugation-closed RHP initialization using unstable plant poles
% and peaks of |Sopt|.

r = floor(m/2);
freqCandidates = [];

if usePu
    pu = pole(Pa);
    pu = pu(real(pu)>0);

    for k=1:numel(pu)
        b = abs(imag(pu(k)));
        if b>0
            freqCandidates(end+1,1)=b; %#ok<AGROW>
        else
            freqCandidates(end+1,1)=max(real(pu(k)),w(1)); %#ok<AGROW>
        end
    end
end

if useS
    crit = abs(Sopt(:));
    pk = false(size(crit));
    if numel(crit)>=3
        pk(2:end-1)=crit(2:end-1)>=crit(1:end-2) & ...
                    crit(2:end-1)>=crit(3:end);
    end
    ids=find(pk);
    [~,ord]=sort(crit(ids),'descend');
    ids=ids(ord);

    for k=1:numel(ids)
        freqCandidates(end+1,1)=w(ids(k)); %#ok<AGROW>
        if numel(freqCandidates)>=3*r
            break
        end
    end
end

% unique separated candidates
freqCandidates = freqCandidates(isfinite(freqCandidates) & freqCandidates>0);
freqCandidates = sort(freqCandidates);

selected=[];
for k=1:numel(freqCandidates)
    if isempty(selected) || ...
       all(abs(log10(freqCandidates(k))-log10(selected))>0.15)
        selected(end+1,1)=freqCandidates(k); %#ok<AGROW>
    end
    if numel(selected)>=r
        break
    end
end

% fill log-spaced if needed
if numel(selected)<r
    fill = logspace(log10(w(1)),log10(w(end)),max(r,2*r)).';
    for k=1:numel(fill)
        if isempty(selected) || ...
           all(abs(log10(fill(k))-log10(selected))>0.10)
            selected(end+1,1)=fill(k); %#ok<AGROW>
        end
        if numel(selected)>=r
            break
        end
    end
end

selected = selected(1:r);

pairs = [];
for k=1:r
    p = realPart + 1i*selected(k);
    pairs(end+1,1)=p; %#ok<AGROW>
    pairs(end+1,1)=conj(p); %#ok<AGROW>
end

if mod(m,2)==1
    spt = [realPart; pairs];
else
    spt = pairs;
end

spt = spt(1:m);
end


function x = pointsToParameters(spt,m,opts)
spt = sanitizeInitialPointSet(spt,m,opts);

if mod(m,2)==0
    r=m/2;
    pos=spt(imag(spt)>0);
    pos=pos(1:r);
    x=[log(real(pos)); log(imag(pos))];
else
    r=(m-1)/2;
    realPts=spt(abs(imag(spt))<1e-12);
    pos=spt(imag(spt)>0);
    x=[log(real(realPts(1))); log(real(pos(1:r))); log(imag(pos(1:r)))];
end
end


function spt = parametersToPoints(x,m,opts)
lbR=opts.PointRealLowerBound;
ubR=opts.PointRealUpperBound;
ubI=opts.PointImagUpperBound;

if mod(m,2)==0
    r=m/2;
    a=exp(x(1:r));
    b=exp(x(r+1:2*r));
    a=min(max(a,lbR),ubR);
    b=min(max(b,lbR),ubI);

    spt=zeros(m,1);
    for k=1:r
        spt(2*k-1)=a(k)+1i*b(k);
        spt(2*k)=a(k)-1i*b(k);
    end
else
    r=(m-1)/2;
    a0=exp(x(1));
    a=exp(x(2:r+1));
    b=exp(x(r+2:2*r+1));

    a0=min(max(a0,lbR),ubR);
    a=min(max(a,lbR),ubR);
    b=min(max(b,lbR),ubI);

    spt=zeros(m,1);
    spt(1)=a0;
    for k=1:r
        spt(2*k)=a(k)+1i*b(k);
        spt(2*k+1)=a(k)-1i*b(k);
    end
end
end


function spt = sanitizeInitialPointSet(spt,m,opts)
spt=spt(:);
if numel(spt)~=m
    error('Initial point count mismatch.');
end

% If user point set is not cleanly conjugation closed, reconstruct a nearby
% conjugate-closed set from representatives.
if mod(m,2)==0
    r=m/2;
    reps=spt(imag(spt)>0);

    if numel(reps)<r
        reps=spt(1:r);
    else
        reps=reps(1:r);
    end

    out=[];
    for k=1:r
        a=max(abs(real(reps(k))),opts.PointRealLowerBound);
        b=max(abs(imag(reps(k))),opts.PointRealLowerBound);
        p=a+1i*b;
        out(end+1,1)=p; %#ok<AGROW>
        out(end+1,1)=conj(p); %#ok<AGROW>
    end
else
    r=(m-1)/2;
    realPts=spt(abs(imag(spt))<1e-12);
    if isempty(realPts)
        a0=max(abs(real(spt(1))),opts.PointRealLowerBound);
    else
        a0=max(abs(real(realPts(1))),opts.PointRealLowerBound);
    end

    reps=spt(imag(spt)>0);
    if numel(reps)<r
        reps=spt(2:min(end,r+1));
    end

    out=a0;
    for k=1:r
        p=reps(min(k,numel(reps)));
        a=max(abs(real(p)),opts.PointRealLowerBound);
        b=max(abs(imag(p)),opts.PointRealLowerBound);
        pp=a+1i*b;
        out(end+1,1)=pp; %#ok<AGROW>
        out(end+1,1)=conj(pp); %#ok<AGROW>
    end
end

spt=out(1:m);
end


function [lb,ub] = parameterBounds(m,opts)
if mod(m,2)==0
    r=m/2;
    lb=[repmat(log(opts.PointRealLowerBound),r,1); ...
        repmat(log(opts.PointRealLowerBound),r,1)];
    ub=[repmat(log(opts.PointRealUpperBound),r,1); ...
        repmat(log(opts.PointImagUpperBound),r,1)];
else
    r=(m-1)/2;
    lb=[log(opts.PointRealLowerBound); ...
        repmat(log(opts.PointRealLowerBound),r,1); ...
        repmat(log(opts.PointRealLowerBound),r,1)];
    ub=[log(opts.PointRealUpperBound); ...
        repmat(log(opts.PointRealUpperBound),r,1); ...
        repmat(log(opts.PointImagUpperBound),r,1)];
end
end


function y = responseOnGrid(G,w)
if isa(G,'frd')
    wr=G.Frequency(:);
    yr=squeeze(G.ResponseData);
    yr=yr(:);

    if numel(wr)==numel(w) && max(abs(wr-w))<1e-12*max(1,max(w))
        y=yr;
    else
        y=interp1(log(wr),yr,log(w),'linear','extrap');
    end
else
    y=squeeze(freqresp(G,w));
    y=y(:);
end
end


function y = evalSymbolic(fun,z)
try
    y=fun(z);
    if numel(y)~=numel(z)
        error('nonvectorized');
    end
catch
    y=arrayfun(fun,z);
end
y=double(y(:));
end


function [Er,Ar,Br,Cr,resid] = realifyGeneral(E,A,B,C,spt)
% Real basis for one real point and/or conjugate pairs.

n=numel(spt);
used=false(n,1);
V=zeros(n,n);
col=0;

for k=1:n
    if used(k), continue, end
    p=spt(k);

    if abs(imag(p))<1e-10
        col=col+1;
        V(k,col)=1;
        used(k)=true;
    else
        [d,j]=min(abs(spt-conj(p)));
        if d>1e-7*max(1,abs(p))
            error('Interpolation set is not conjugation closed.');
        end

        col=col+1;
        V(k,col)=1/sqrt(2);
        V(j,col)=1/sqrt(2);

        col=col+1;
        V(k,col)=-1i/sqrt(2);
        V(j,col)= 1i/sqrt(2);

        used(k)=true;
        used(j)=true;
    end
end

Ec=V'*E*V;
Ac=V'*A*V;
Bc=V'*B;
Cc=C*V;

scale=max(1,max([norm(Ec,'fro'),norm(Ac,'fro'),norm(Bc),norm(Cc)]));
resid=max([norm(imag(Ec),'fro'),norm(imag(Ac),'fro'), ...
           norm(imag(Bc)),norm(imag(Cc))])/scale;

Er=real(Ec);
Ar=real(Ac);
Br=real(Bc);
Cr=real(Cc);
end
