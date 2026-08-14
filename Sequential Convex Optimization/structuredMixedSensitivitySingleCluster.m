function [C,info] = structuredMixedSensitivitySingleCluster(P,Copt,W1,W3,m,relativeOrder,X,opts)
% STRUCTUREDMIXEDSENSITIVITYSINGLECLUSTER
%
%   [C,info] = structuredMixedSensitivitySingleCluster( ...
%       P,Copt,W1,W3,m,relativeOrder,X,opts)
%
% Reduced-order mixed-sensitivity controller approximation with NO supplied
% reduced-order initial controller.
%
% Inputs
%   P             SISO FRD (or struct Frequency/ResponseData), performance plant
%   Copt          SISO FRD/LTI optimal high-order controller response
%   W1,W3         sensitivity / complementary-sensitivity weights
%   m             desired denominator order of reduced controller C
%   relativeOrder imposed denominator-minus-numerator relative degree
%   X             rational SISO plant approximation used ONLY for allmargin
%   opts          options; opts.r is the degree of T(s)=(s-z)^r
%
% True sampled merit:
%   gamma(C)=max_w sqrt(|W1*S|^2+|W3*P*C*S|^2), S=1/(1+P*C).
%
% Reference-oriented approximation residual:
%   Eapp = Wx*P*(Copt-C)*S*Sopt,
%   Sopt = 1/(1+P*Copt),  Wx = sqrt(|W1|^2+|W3|^2).
%
% Initialization uses opts.C0 when supplied; otherwise weightedSsest fits the
% FRD Copt with identification weight |Wx*P*Sopt^2|.  In either case the
% initializer denominator is retained and the requested numerator degree is
% enforced before the main sequential SOCP.  The SOCP is softly guided at
% active approximation-error peaks by T(s)=(s-z)^r.
%
% The SOCP uses at most opts.maxOptSamples points; true gamma and line-search
% acceptance always use the complete FRD grid.
%
% Requires YALMIP + SDPT3.

if nargin < 6 || isempty(relativeOrder), relativeOrder = 0; end
if nargin < 7 || isempty(X)
    error('X must be supplied as a rational SISO plant approximation.');
end
if nargin < 8 || isempty(opts), opts = struct; end
opts = defaultOptions(opts);

validateattributes(m,{'numeric'},{'scalar','integer','positive'});
validateattributes(relativeOrder,{'numeric'},{'scalar','integer','nonnegative'});
if relativeOrder > m
    error('relativeOrder cannot exceed m.');
end
numDegree = m-relativeOrder;

assert(issiso(X),'X must be SISO.');
if isa(X,'frd'), error('X must be rational tf/zpk/ss, not FRD.'); end
X = tf(X);

[w,Pj] = unpackFRD(P);
[w,ord] = sort(w(:)); Pj = Pj(ord);
[w,ia] = unique(w,'stable'); Pj = Pj(ia);
if any(w<0) || numel(w)<4, error('P requires >=4 nonnegative frequency samples.'); end

Coptj = evaluateResponse(Copt,w);
W1j = evaluateWeight(W1,w);
W3j = evaluateWeight(W3,w);
Sopt = 1./(1+Pj.*Coptj);
Wxj = sqrt(abs(W1j).^2 + abs(W3j).^2);

if any(~isfinite(Sopt))
    error('Copt produces nonfinite Sopt on the supplied P frequency grid.');
end

r = opts.r;
if isempty(r)
    error(['For FRD Copt, opts.r must be supplied explicitly for ', ...
           'T(s)=(s-z)^r (e.g. opts.r=n-m-1).']);
end
validateattributes(r,{'numeric'},{'scalar','integer','nonnegative'});

% z search interval is frequency-grid based because Copt may be pure FRD.
wp = w(w>0);
if isempty(wp), wp = 1; end
zmin = max(opts.absoluteZmin, opts.zMinFactor*min(wp));
zmax = max(zmin*(1+1e-8), opts.zMaxFactor*max(wp));
if ~isempty(opts.zRange)
    zmin = opts.zRange(1); zmax = opts.zRange(2);
end

fprintf('\nSingle-cluster mixed-sensitivity reduction from Copt FRD\n');
fprintf('========================================================\n');
fprintf('reduced denominator order = %d\n',m);
fprintf('relative order            = %d\n',relativeOrder);
fprintf('numerator degree <=       = %d\n',numDegree);
fprintf('T degree r                = %d\n',r);
fprintf('frequency samples         = %d\n',numel(w));
fprintf('z range                    = [%.4e, %.4e]\n',zmin,zmax);

%% Initialization: opts.C0 if supplied, otherwise weightedSsest
if ~isempty(opts.C0)
    if isa(opts.C0,'frd') || ~issiso(opts.C0)
        error('opts.C0 must be a rational SISO tf/zpk/ss model.');
    end

    Craw = ss(opts.C0);
    initMethod = 'opts.C0';
else
    Winitj = abs(Wxj.*Pj.*Sopt.^2);
    if max(Winitj) <= 0 || any(~isfinite(Winitj))
        Winitj = ones(size(Winitj));
    end

    CoptFRD = frd(Coptj,w);
    Winit = frd(Winitj,w);
    Craw = weightedSsest(CoptFRD,Winit,m);
    initMethod = 'weightedSsest';
end

[n0,d0] = tfdata(tf(Craw),'v');
n0 = real(n0(:).');
d0 = real(d0(:).');
a = d0(1); n0 = n0/a; d0 = d0/a;

m0 = length(d0)-1;
if m0 ~= m
    error('%s initializer has denominator order %d, but requested m=%d.', ...
          initMethod,m0,m);
end

% Enforce the requested relative degree while retaining the initializer
% denominator.
n0 = restrictNumeratorDegree(n0,numDegree);
C = tf(n0,d0);

marginInit = allmarginStabilityTest(X,C);
[gammaInit,~,~,stateInit] = mixedCostSamples(Pj,W1j,W3j,w,C,opts);
approxInit = approximationCost(Pj,Coptj,Wxj,Sopt,w,C);

if ~marginInit.stable || ~stateInit.valid
    error(['%s initializer does not provide an admissible stabilizing ', ...
           'order-%d controller. allmargin reason: %s'], ...
          initMethod,m,marginInit.reason);
end

fprintf('initializer                = %s\n',initMethod);
fprintf('initial true gamma         = %.12e\n',gammaInit);
fprintf('initial approximation norm = %.12e\n\n',approxInit);

%% Main sequential SOCP
Niter = opts.maxIter;
info.gamma = nan(Niter+1,1);
info.approxGamma = nan(Niter+1,1);
info.z = nan(Niter,1);
info.TfitResidual = nan(Niter,1);
info.alpha = nan(Niter,1);
info.dGamma = nan(Niter,1);
info.tStructure = nan(Niter,1);
info.activeW = cell(Niter,1);
info.initialization.method = initMethod;
info.initialization.gamma = gammaInit;
info.initialization.approxGamma = approxInit;
info.initialization.margins = marginInit;
info.Sopt = Sopt;
info.Wx = Wxj;
info.r = r;

[gamma,~,gVec,state] = mixedCostSamples(Pj,W1j,W3j,w,C,opts);
info.gamma(1)=gamma;
info.approxGamma(1)=approximationCost(Pj,Coptj,Wxj,Sopt,w,C);
trustN=opts.trustN; trustD=opts.trustD;

for iter=1:Niter
    [nC,dC]=tfdata(tf(C),'v');
    nC=real(nC(:).'); dC=real(dC(:).');
    a=dC(1); nC=nC/a; dC=dC/a;
    nC=restrictNumeratorDegree(nC,numDegree);

    [gamma,wWorst,gVec,state]=mixedCostSamples(Pj,W1j,W3j,w,C,opts);
    Cj=evalLTI(C,w);
    Sa=1./(1+Pj.*Cj);
    Eapp=Wxj.*Pj.*(Coptj-Cj).*Sa.*Sopt;
    appMag=abs(Eapp);

    nActive=m+1;
    if ~isempty(opts.numActive), nActive=opts.numActive; end
    [wA,idxA]=activePeaksSamples(w,appMag,nActive,opts);
    info.activeW{iter}=wA;

    % Frozen-weight approximation gradient wrt C: dEapp ~= -Wapp*dC.
    Wapp=Wxj.*Pj.*Sa.*Sopt;
    Happ=-conj(Eapp).*Wapp./max(abs(Eapp),1e-14);
    Tfit=fitSingleClusterTFRD(Happ(idxA),dC,wA,r,zmin,zmax,opts);

    info.z(iter)=Tfit.z;
    info.TfitResidual(iter)=Tfit.residual;

    fprintf('Iteration %2d: gamma %.12e, app %.12e, z %.5e, Tres %.4e\n', ...
        iter,gamma,max(appMag),Tfit.z,Tfit.residual);

    % Build a reduced SOCP grid, always retaining all active points.
    idxOpt = optimizationGridIndices(w,idxA,opts.maxOptSamples);
    [~,idxAopt] = ismember(idxA,idxOpt);

    [dn,da,dGamma,step]=mixedSensitivityDirection( ...
        Pj(idxOpt),W1j(idxOpt),W3j(idxOpt),w(idxOpt), ...
        nC,dC,numDegree,gamma,idxAopt,Tfit,trustN,trustD,opts);
    info.dGamma(iter)=dGamma;
    if isfield(step,'tStructure'), info.tStructure(iter)=step.tStructure; end

    if ~step.success
        trustN=trustN*opts.trustShrink; trustD=trustD*opts.trustShrink;
        info.gamma(iter+1)=gamma; info.approxGamma(iter+1)=max(appMag);
        if trustD<opts.minTrust, break; end
        continue
    end

    accepted=false; alpha=1;
    while alpha>=opts.alphaMin
        nt=nC+alpha*dn(:).';
        dt=dC; dt(2:end)=dt(2:end)+alpha*da(:).';
        if any(real(roots(dt))>=-opts.stabilityMargin)
            alpha=alpha/2; continue
        end
        Ct=tf(nt,dt);
        [gt,~,~,stt]=mixedCostSamples(Pj,W1j,W3j,w,Ct,opts);
        mt=allmarginStabilityTest(X,Ct);
        if mt.stable && stt.valid && gt < gamma-opts.improvementTol*max(1,abs(gamma))
            accepted=true; break
        end
        alpha=alpha/2;
    end
    info.alpha(iter)=alpha;

    if accepted
        old=gamma; C=Ct; gamma=gt;
        info.gamma(iter+1)=gamma;
        info.approxGamma(iter+1)=approximationCost(Pj,Coptj,Wxj,Sopt,w,C);
        trustN=min(opts.maxTrust,trustN*opts.trustGrow);
        trustD=min(opts.maxTrust,trustD*opts.trustGrow);
        if (old-gamma)/max(old,eps)<opts.relTol, break; end
    else
        info.gamma(iter+1)=gamma; info.approxGamma(iter+1)=max(appMag);
        trustN=trustN*opts.trustShrink; trustD=trustD*opts.trustShrink;
        if trustD<opts.minTrust, break; end
    end
end

[gamma,wWorst,~,state]=mixedCostSamples(Pj,W1j,W3j,w,C,opts);
k=find(~isnan(info.gamma),1,'last'); if isempty(k), k=1; end
info.gamma=info.gamma(1:k);
info.approxGamma=info.approxGamma(1:k);
info.iterations=k-1;
info.finalGamma=gamma;
info.finalApproxGamma=approximationCost(Pj,Coptj,Wxj,Sopt,w,C);
info.finalC=C;
info.finalWorstFrequency=wWorst;
info.finalMinReturnDifference=state.minReturnDifference;
info.finalMargins=allmarginStabilityTest(X,C);
info.frequency=w;
info.Copt=Copt;
info.X=X;
info.m=m;
info.relativeOrder=relativeOrder;
end

function val=approximationCost(Pj,Coptj,Wxj,Sopt,w,C)
Cj=evalLTI(C,w); Sa=1./(1+Pj.*Cj);
val=max(abs(Wxj.*Pj.*(Coptj-Cj).*Sa.*Sopt));
end

function Tout=fitSingleClusterTFRD(HA,dC,wA,r,zmin,zmax,opts)
s=1j*wA(:); HA=HA(:);
% With Copt available only as FRD, no high-order denominator polynomial is
% assumed. The current reduced denominator contributes the known phase.
A=polyval(dC,-s)./polyval(dC,s);
R=normalizePhasor(HA)./normalizePhasor(A); R=normalizePhasor(R);
weights=abs(HA);
if max(weights)<=0 || ~all(isfinite(weights)), weights=ones(size(weights));
else, weights=weights/max(weights); end
weights=max(weights.^opts.TfitWeightPower,1e-8);
if r==0
    U=ones(size(s)); z=sqrt(zmin*zmax);
else
    fun=@(lz) singleResidual(exp(min(max(lz,log(zmin)),log(zmax))),R,s,weights,r);
    lgrid=linspace(log(zmin),log(zmax),opts.Tgrid);
    vals=arrayfun(fun,lgrid); [~,ii]=min(vals);
    lz=fminsearch(fun,lgrid(ii),optimset('Display','off','TolX',opts.TfitTol, ...
        'TolFun',opts.TfitTol,'MaxIter',opts.TfitMaxIter));
    z=exp(min(max(lz,log(zmin)),log(zmax)));
    U=((s+z)./(-s+z)).^r;
end
Tout.z=z; Tout.residual=weightedTResidual(R,U,weights);
Tout.targetH=normalizePhasor(A.*U);
Tout.T=poly(z*ones(1,r));
end

function v=singleResidual(z,R,s,wgt,r)
U=((s+z)./(-s+z)).^r; v=weightedTResidual(R,U,wgt);
end
function val=weightedTResidual(R,U,wgt)
val=sum(wgt.*abs(R-U).^2)/sum(wgt);
end

function [dn,da,dGamma,out]=mixedSensitivityDirection( ...
    Pj,W1j,W3j,w,nC,dC,numDegree,gamma,idxA,Tfit,trustN,trustD,opts)
m=length(dC)-1; out.success=false;
dn=zeros(numDegree+1,1); da=zeros(m,1); dGamma=0;
scaleN=max(abs(nC(:)),1); scaleD=max(abs(dC(2:end).'),1);
uN=sdpvar(numDegree+1,1); uD=sdpvar(m,1); dg=sdpvar(1); ts=sdpvar(1);
s=1j*w(:); dv=polyval(dC,s); nv=polyval(nC,s); Cj=nv./dv;
S=1./(1+Pj.*Cj); F1=W1j.*S; F3=W3j.*Pj.*Cj.*S;
Bn=polynomialBasis(s,numDegree); Bd=polynomialBasis(s,m-1);
JN=(Bn./dv).*scaleN.'; JD=(-(nv./dv.^2).*Bd).*scaleD.';
G1=-W1j.*Pj.*S.^2; G3=W3j.*Pj.*S.^2;
Cons=[gamma+dg>=0, dg<=opts.maxPositivePredictedGamma*gamma, ...
      dg>=-opts.maxRelativeGammaStep*gamma];
for k=1:numel(w)
    dqN=JN(k,:); dqD=JD(k,:);
    f1=F1(k)+G1(k)*(dqN*uN+dqD*uD);
    f3=F3(k)+G3(k)*(dqN*uN+dqD*uD);
    Cons = [Cons, cone([real(f1); imag(f1); real(f3); imag(f3)],gamma+dg)]; %#ok<AGROW>
end
Cons = [Cons, cone(uN,trustN), cone(uD,trustD), ts>=0];
for ii=1:numel(idxA)
    k=idxA(ii); dq=JN(k,:)*uN+JD(k,:)*uD;
    dir=-conj(Tfit.targetH(ii)); rot=conj(dir)*dq;
    sc=max(abs(Cj(k)),1e-6);
    yi = imag(rot)/sc;
    Cons = [Cons, yi <= ts, -yi <= ts, ...
        real(rot)/sc >= -opts.structureDirectionSlack*ts]; %#ok<AGROW>
end
obj=dg/max(gamma,1e-12)+opts.extremalWeight*ts;
sol=optimize(Cons,obj,sdpsettings('solver','sdpt3','verbose',opts.solverVerbose));
if sol.problem~=0, out.problem=sol.problem; out.info=sol.info; return; end
un=value(uN); ud=value(uD);
if any(~isfinite(un)) || any(~isfinite(ud)), return; end
dn=(scaleN.*un).'; dn=dn(:); da=scaleD.*ud;
dGamma=value(dg); out.success=true; out.tStructure=value(ts);
end

function [gamma,wWorst,g,state]=mixedCostSamples(Pj,W1j,W3j,w,C,opts)
Cj=evalLTI(C,w); D=1+Pj.*Cj; S=1./D;
F1=W1j.*S; F3=W3j.*Pj.*Cj.*S;
g=sqrt(abs(F1).^2+abs(F3).^2);
if any(~isfinite(g)) || any(~isfinite(D)), gamma=inf; wWorst=NaN;
else, [gamma,ii]=max(g); wWorst=w(ii); end
state.minReturnDifference=min(abs(D));
state.valid=isfinite(gamma) && state.minReturnDifference>=opts.returnDifferenceMargin;
end

function idx=optimizationGridIndices(w,idxActive,maxN)
% Reduced grid for the convex SOCP. Full FRD grid is still used for
% performance evaluation and line-search acceptance.
N=numel(w);
idxActive=unique(idxActive(:));

if isempty(maxN) || maxN<=0 || N<=maxN
    idx=(1:N).';
    return
end

maxN=max(round(maxN),numel(idxActive)+2);

% Index-uniform sampling works well when the supplied FRD grid is already
% logarithmic (the usual case). Active points are added explicitly.
base=unique(round(linspace(1,N,maxN-numel(idxActive))));
idx=unique([base(:); idxActive; 1; N]);

% If adding active points pushed us above maxN, retain all active points and
% thin only the non-active background samples.
if numel(idx)>maxN
    nonActive=setdiff(idx,idxActive,'stable');
    keep=maxN-numel(idxActive);
    if keep>0
        jj=unique(round(linspace(1,numel(nonActive),keep)));
        nonActive=nonActive(jj);
    else
        nonActive=[];
    end
    idx=unique([idxActive; nonActive(:)]);
end
idx=sort(idx);
end

function [wA,idxA]=activePeaksSamples(w,g,nWanted,opts)
N=numel(w); idx=[];
if N>=3, idx=find(g(2:end-1)>=g(1:end-2)&g(2:end-1)>=g(3:end))+1; end
if N>=2 && g(1)>=g(2), idx=[1;idx(:)]; end
if N>=2 && g(end)>=g(end-1), idx=[idx(:);N]; end
if isempty(idx), [~,ii]=max(g); idx=ii; end
mp=g(idx); gm=max(mp); ia=idx(mp>=opts.peakThreshold*gm);
if numel(ia)<nWanted
    [~,o]=sort(mp,'descend'); ia=idx(o(1:min(nWanted,numel(o))));
end
if numel(ia)>nWanted
    [~,o]=sort(g(ia),'descend'); ia=ia(o(1:nWanted));
end
idxA=sort(ia(:)); wA=w(idxA);
end

function [w,H]=unpackFRD(G)
if isa(G,'frd')
    w=G.Frequency(:); H=reshape(squeeze(G.ResponseData),[],1);
elseif isstruct(G)&&isfield(G,'Frequency')&&isfield(G,'ResponseData')
    w=G.Frequency(:); H=reshape(squeeze(G.ResponseData),[],1);
else
    error('P must be SISO FRD or Frequency/ResponseData struct.');
end
if numel(w)~=numel(H), error('P frequency/response sizes disagree.'); end
end

function H=evaluateResponse(G,w)
if isa(G,'frd') || isa(G,'tf') || isa(G,'zpk') || isa(G,'ss')
    H=evalLTI(G,w);
elseif isstruct(G)&&isfield(G,'Frequency')&&isfield(G,'ResponseData')
    wf=G.Frequency(:); hf=reshape(squeeze(G.ResponseData),[],1);
    H=interp1(wf,hf,w,'linear','extrap');
else
    error('Copt must be SISO FRD/LTI or Frequency/ResponseData struct.');
end
H=H(:);
end

function H=evaluateWeight(W,w)
if isnumeric(W)&&isscalar(W), H=repmat(W,numel(w),1);
else, H=evaluateResponse(W,w); end
end

function H=evalLTI(G,w)
R=freqresp(G,w(:)); H=reshape(R,[],1);
end

function out=allmarginStabilityTest(X,C)
out.stable=false; out.margins=[]; out.reason='';
try
    M=allmargin(X*C); out.margins=M;
    if isempty(M)||~isfield(M,'Stable'), out.reason='No Stable field'; return; end
    f=[M.Stable]; out.stable=all(isfinite(f))&&all(logical(f));
    if out.stable, out.reason='stable according to allmargin(X*C)';
    else, out.reason='allmargin(X*C).Stable is false'; end
catch ME
    out.reason=['allmargin failed: ' ME.message];
end
end

function n=restrictNumeratorDegree(n,numDegree)
n=real(n(:).'); N=numDegree+1;
if numel(n)>N, n=n(end-N+1:end);
elseif numel(n)<N, n=[zeros(1,N-numel(n)),n]; end
end

function B=polynomialBasis(s,degree)
s=s(:); B=zeros(numel(s),degree+1);
for k=0:degree, B(:,k+1)=s.^(degree-k); end
end

function y=normalizePhasor(x)
x=x(:); y=ones(size(x)); a=abs(x);
ii=isfinite(real(x))&isfinite(imag(x))&(a>1e-14); y(ii)=x(ii)./a(ii);
end

function o=defaultOptions(o)
o=setdef(o,'C0',[]);
o=setdef(o,'r',[]);
o=setdef(o,'maxIter',30); o=setdef(o,'maxOptSamples',4000); o=setdef(o,'relTol',1e-6);
o=setdef(o,'improvementTol',1e-8); o=setdef(o,'numActive',[]);
o=setdef(o,'peakThreshold',0.97);
o=setdef(o,'Tgrid',1e3); o=setdef(o,'TfitTol',1e-7);
o=setdef(o,'TfitMaxIter',200); o=setdef(o,'TfitWeightPower',2);
o=setdef(o,'absoluteZmin',1e-6); o=setdef(o,'zMinFactor',1);
o=setdef(o,'zMaxFactor',10); o=setdef(o,'zRange',[]);
o=setdef(o,'extremalWeight',0.001); o=setdef(o,'structureDirectionSlack',0.05);
o=setdef(o,'trustN',0.20); o=setdef(o,'trustD',0.10);
o=setdef(o,'minTrust',1e-4); o=setdef(o,'maxTrust',1);
o=setdef(o,'trustGrow',1.5); o=setdef(o,'trustShrink',0.5);
o=setdef(o,'maxRelativeGammaStep',0.1);
o=setdef(o,'maxPositivePredictedGamma',0.1);
o=setdef(o,'alphaMin',1/128); o=setdef(o,'stabilityMargin',1e-8);
o=setdef(o,'returnDifferenceMargin',1e-5); o=setdef(o,'solverVerbose',0);
end
function s=setdef(s,name,val)
if ~isfield(s,name)||isempty(s.(name)), s.(name)=val; end
end
