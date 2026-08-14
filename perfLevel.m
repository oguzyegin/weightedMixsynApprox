function gamma = perfLevel(P,C,W1,W2,w)


S= 1/(1+C*P);
T= C*P/(1+C*P);

e1 = squeeze(freqresp(W1*S,w));
e2 = squeeze(freqresp(W2*T,w));

gamma = sqrt(abs(e1).^2+abs(e2).^2);
