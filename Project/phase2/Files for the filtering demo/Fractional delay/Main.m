% Fractional delay of sampled data
% Obtain the samples of a delayed band-limited continuous time signal.

% Oktay Sipahigil
% 02/01/2017

%% Samples of a bandlimited signal
xfcn = @(t)cos(2*pi*2*t) + 0.4*cos(2*pi*3*t+0.8);
T = 0.15;
dur = 10;
fsplot = 1000;
t = (0:fsplot*dur-1)/fsplot;
n = 0:dur/T-1;
x = xfcn(n*T);
figure;
plot(t,xfcn(t),'--'); hold on; stem(n*T,x,'filled'); xlabel('Time (sec)');
legend('CT signal','Its sampled version');

%% Fractional delay filter
td = 0.05;
nd = td/T;
N = 41;
nn = -(N-1)/2 : (N-1)/2;
win = hamming(N)';
h = win.*sinc(nn-nd);
figure; stem(nn,h,'filled'); xlabel('n'); ylabel('h[n]');

%% Result
y = filter(h,1,x);
figure;
ha1 = subplot(2,1,1); plot(t,xfcn(t),'--'); hold on; stem(n*T,x,'filled'); xlabel('Time (sec)');
title('x(t) and samples');
ha2 = subplot(2,1,2); plot(t,xfcn(t-td),'--'); hold on; stem(n((N+1)/2:end)*T,y((N+1)/2:end),'filled'); xlabel('Time (sec)');
title('x(t-t_d) and samples');
linkaxes([ha1 ha2]);