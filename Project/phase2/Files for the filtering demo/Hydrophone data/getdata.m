function [x,fs] = getdata()
seed = 1;
rngs = rng(seed);   onCleanupHandle = onCleanup(@()rng(rngs));

fs = 10e3;
dur = 3;

t = (0:dur*fs-1)/fs;
x = 0.01*randn(size(t));
f = 0.5;   a = 1;     x = x + a*cos(2*pi*f*t + 2*pi*rand(1));
f = 1.9;   a = 0.2;   x = x + a*cos(2*pi*f*t + 2*pi*rand(1));

t0 = 1;   d = 0.1;   i = t0<=t&t<t0+d;
f = 0.5e3;   a = 0.2;   x(i) = x(i) + a*cos(2*pi*f*t(i) + 2*pi*rand(1));

t0 = 1.4;   d = 0.05;   i = t0<=t&t<t0+d;
f = 2e3;   a = 0.2;   x(i) = x(i) + a*cos(2*pi*f*t(i) + 2*pi*rand(1));

t0 = 2.2;   d = 0.15;   i = t0<=t&t<t0+d;
f = 0.7e3;   a = 0.2;   x(i) = x(i) + a*cos(2*pi*f*t(i) + 2*pi*rand(1));

x = x + 0.1*cos(2*pi*0.95e3*t + 2*pi*rand(1));

x = 10*x;