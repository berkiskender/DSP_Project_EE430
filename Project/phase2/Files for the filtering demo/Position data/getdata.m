function [x,fs] = getdata()
seed = 6;
rngs = rng(seed);   onCleanupHandle = onCleanup(@()rng(rngs));

fs = 100;
dur = 3;
t = (0:dur*fs-1)/fs;
x = 0.01*randn(size(t));

i = 0<=t&t<1;   x(i) = x(i) + t(i).^2;
i = 1<=t&t<2;   x(i) = x(i) + 0.1*t(i) + 0.9;
i = 2<=t&t<3;   x(i) = x(i) - 0.3*t(i) + 1.7;
