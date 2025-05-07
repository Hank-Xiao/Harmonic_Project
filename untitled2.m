%实现对原始信号的插值
% rng default
% nominalFs = 1024;
% f = 50;
% Tx = 0:0.02/(nominalFs-1):0.02;
% irregTx = sort(Tx + 1e-4*rand(size(Tx)));
% x = 100*cos(2*pi*f*Tx)+20*cos(2*pi*f*2*Tx)+10*cos(2*pi*f*3*Tx)+5*cos(2*pi*f*4*Tx)+cos(2*pi*f*5*Tx);
% plot(Tx,x,'.')


