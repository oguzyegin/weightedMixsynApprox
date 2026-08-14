function [Cf,gamma0,gamma] = callSystune(X,W1,W2,Qi,wHF)

if nargin < 5 || isempty(wHF)
    wHF = 1e6;
end

s = tf('s');

X  = tf(X);
W1 = tf(W1);
W2 = tf(W2);
Qi = minreal(tf(Qi));

%% Properize W2 if necessary

q = max(0,(length(find(W2.num{1}))-1)-(length(find(W2.den{1}))-1));

W2app = W2*(wHF/(s+wHF))^q;

%% Tunable controller initialized exactly with Qi
Qi = tf(Qi);
nQi = Qi.num{1}; nK = nQi(find(nQi));
dQi = Qi.den{1}; dK = dQi(find(dQi));

K = tunableTF('K',tf(nK,dK));

X.InputName  = 'u';
X.OutputName = 'y';

K.InputName  = 'e';
K.OutputName = 'u';

W1s = ss(W1);
W1s.InputName  = 'e';
W1s.OutputName = 'z1';

W2s = ss(W2app);
W2s.InputName  = 'y';
W2s.OutputName = 'z2';

Sum = sumblk('e = r-y');

%% Generalized closed loop

CL0 = connect(X,K,W1s,W2s,Sum, ...
              'r',{'z1','z2'});

%% Tune

opt1 = systuneOptions('SoftTol',1e-6);
ReqPerf = TuningGoal.Gain('r',{'z1','z2'},1);
% Force controller K to be stable
ReqKstable = TuningGoal.ControllerPoles('K',0,0,Inf);

opt1.Hidden.Problem = 'Hinf';
opt1.Display = 'final';
opt1.MaxIter = 1e3;
CL = systune(CL0,ReqPerf,ReqKstable,opt1);

Cf = tf(CL.Blocks.K);


end