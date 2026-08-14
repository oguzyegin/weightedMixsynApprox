function X = weightedSsest(P,W,n)
% WEIGHTEDSSEST  Weighted order-n state-space fit of FRD P.
%
% This helper is the user-specified initializer.

w = P.Frequency(:);
H = squeeze(freqresp(P,w));

data = idfrd(H,w,0);

Wk = abs(squeeze(freqresp(W,w)));
Wk = Wk(:)/max(Wk);

opt = ssestOptions;
opt.WeightingFilter = Wk;
opt.Display = 'off';

X = ssest(data,n,opt);
X = ss(X);
end
