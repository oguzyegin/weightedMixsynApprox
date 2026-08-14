clear all
close all
clc

%%
% Define weights and the plant
s = tf('s');
W1 = 0.5/(s+1e-5);
W2 = (0.01*s^2+0.3*s+0.2);

nP = exp(-0.2*s) * (5-s*exp(-0.1*s));
dP = ((s+1)*(3*s+0.5)+(2*s+7)*exp(-1.5*s)+(s-1)*exp(-2*s));
P = nP/dP;
%% 
addpath(genpath([pwd, '/Hinf Optimal with YALTA v3']))
% Convert plant parameters into YALTA input
[N,D,tau_n,tau_d,T0,~] = nPdP2YALTA(nP,dP);

% Detect unstable poles by using YALTA
scale = 1e9;
v = round(tau_d*scale);
g = v(1);
for k = 2:numel(v), g = gcd(g,v(k)); end
iTau = g/scale;
iDelayVect = round(tau_d/iTau);
Z = delayFrequencyAnalysisMin(D,iDelayVect(2:end)',1,iTau);
nMd = real(poly(Z.UnstablePoles)); dMd = real(poly(-Z.UnstablePoles));
Md = tf(nMd,dMd);
disp('Unstable poles of the given plant:');
Z.UnstablePoles

%% 
% Complete rest of the decomposition by using symbolic variables
% P = Mn*No/Md
close all;
syms s;
q = ((s+1)*(3*s+0.5)+(2*s+7)*exp(-1.5*s)+(s-1)*exp(-2*s));

Mns = exp(-T0*s)* (5-s*exp(-0.1*s))/ (s+5*exp(-0.1*s));
Mds = poly2sym(nMd,s)/poly2sym(dMd,s);
Nos = (s+5*exp(-0.1*s)) *Mds/q;

%% 
% Find inf. dim. optimal controller

nw1 = W1.num{1};
dw1 = W1.den{1};
nw2 = W2.num{1};
dw2 = W2.den{1};
warning off;
[Copt, g, ~, ~, ~] = TO_hinfsyn(nw1,dw1,nw2,dw2,Z.UnstablePoles,Mns,Mds,Nos);
warning on;

%% 
% Frequency response data of P and Copt
close all;
warning off;
w = logspace(-5,3,1e3);
s = 1j*w;
Pfr = squeeze(freqresp(P,w));
Cfr = eval(Copt).';
Pw = frd(Pfr,w);
Cw = frd(Cfr,w);
warning on;

%%
% Find an m-th order inital Controller by using ssest 
m = 8; % desired order for approximation
W1w = squeeze(freqresp(W1,w));
W2w = squeeze(freqresp(W2,w));
Sopt = 1./(1+Pfr.*Cfr);

Wx = sqrt(abs(W1w).^2 + abs(W2w).^2);
Winit = abs(Wx.*Pfr.*Sopt.^2);
Winit = frd(Winit,w);
addpath(fullfile(pwd, 'Sequential Convex Optimization'))
tic
Ca = weightedSsest(Cw,Winit,m);
t_ss = toc;
Ca_30 = weightedSsest(Cw,Winit,30); % for hifoo and irka

%%
% Controller designed by Vector fitting
nC = [51.988, 1.1223e6, 1.5798e7, 6.7827e7, 3.3433e8, 3.6536e8, 8.8553e7, 1.8972e8];
dC = [1, 6365.80, 1.8877e5, 2.2990e6, 1.7066e7, 8.5950e7, 2.2647e8, 2.9232e8, 554.24];
% tic
% [n0, d0, ~, ~] = vf_fit_order8(w, Cfr, 8);
% [nC, dC, ~, ~] = vf_tune_hinf(n0, d0, w);
Cvf = tf(nC,dC);
% t_vf = toc;

%% Call systune
% Find a high order approximation of the plant, Pa, to be used in systune

opt = tfestOptions('EnforceStability',true);

Pa = tfest(frd(squeeze(freqresp(P*Md,w)),w),15,14,opt);
Pa = minreal(ss(Pa/Md),1e-6);

% Systune
tic
Cs = callSystune(Pa,W1,W2,Ca);
t_s = toc;
% Systune with vf initialization
tic
Csv = callSystune(Pa,W1,W2,Cvf);
t_sv = toc;

%% Call sequential convex minimization by proposed T structure
addpath(genpath(fullfile(pwd,'YALMIP')));
yalmip('clear')
relOrd = 1;
opts = struct;
opts.r = 2*m+2; % Order of T for inf dim Copt
opts.C0 = Ca; % Weighted ssest returns initial controller
tic
[Cp,info_p] = structuredMixedSensitivitySingleCluster(Pw,Cw,W1,W2,m,relOrd,Pa,opts);
t_p = toc;
opts.C0 = Cvf; % Vf returns initial controller
tic
[Cpv,info_pv] = structuredMixedSensitivitySingleCluster(Pw,Cw,W1,W2,m,relOrd,Pa,opts);
t_pv = toc;

%% TFIRKA based mixed s min
% Initial interpolation points: unst poles + 
% nyquist encirclement intersections on real axis
% pu = pole(Pa);
%     pu = pu(real(pu) > 0);
% 
%     opts.InitialPoints = [ ...
%         pu.', ...
%         0.02+0.7j, 0.02-0.7j, ...
%         0.1-4.2j, 0.1+4.2j, ...
%         1-15j,   1+15j ];
tic
[Ci,~,info_irka] = tfirkaMixedSensitivity(Copt,Pw,Pa,W1,W2,m);
t_i = toc;
%% TDS-control
warning off;
addpath(genpath([pwd, '\tds-control-main']));
tic
[Ct,info_tds] = designTDSMixedSensitivity8(Ca);
t_t = toc;
tic
[Ctv,info_tdsv] = designTDSMixedSensitivity8(Cvf);
t_tv = toc;
warning on;

%% HIFOO
% addpath(genpath([pwd, '\HIFOO3.501']));
% addpath(genpath([pwd, '\hanso2_02']));
% [Ch,gamma,info_h] = hifooMixedSensitivityMOR(Ca_30,Pw,Pa,ss(W1),ss(W2p),m,Ca);
% [Chv,~,~] = hifooMixedSensitivityMOR(Ca_30,Pw,Pa,ss(W1),ss(W2),m,Cvf);

%%
% Performance Levels
wp = logspace(-5,4.3,1e5);
gss = perfLevel(P,Ca,W1,W2,wp);
gvf = perfLevel(P,Cvf,W1,W2,wp);
gs = perfLevel(P,Cs,W1,W2,wp);
gsv = perfLevel(P,Csv,W1,W2,wp);
gt = perfLevel(P,Ct,W1,W2,wp);
gtv = perfLevel(P,Ctv,W1,W2,wp);
gi = perfLevel(P,Ci,W1,W2,wp);
gp = perfLevel(P,Cp,W1,W2,wp);
gpv = perfLevel(P,Cpv,W1,W2,wp);

Method = [
    "Ssest"
    "Vector fitting"
    "TF-irka"
    "systune (ssest init.)"
    "systune (VF init.)"
    "TDS-control  (ssest init.)"
    "TDS-control (VF init.)"
    "Proposed (ssest init.)"
    "Proposed (VF init.)"
];

MaxPerformance = [
    max(abs(gss))
    max(abs(gvf))
    max(abs(gi))
    max(abs(gs))
    max(abs(gsv))
    max(abs(gt))
    max(abs(gtv))
    max(abs(gp))
    max(abs(gpv))
];
TimeCost = [
    t_ss
    t_vf
    t_i
    t_s
    t_sv
    t_t
    t_tv
    t_p
    t_pv
];

Results = table(Method, MaxPerformance, TimeCost);

disp(Results)

figure(1); clf;
semilogx(wp,gtv); hold on;
semilogx(wp,gpv);