function [x,fs,t] = getAudio(file,t0,t1)

info = audioinfo(file);
fs = info.SampleRate;
N  = info.TotalSamples;

n0 = t0*fs;
n1 = t1*fs-1;
n1 = min(n1,N-1);

x = audioread(file,[n0 n1]+1);   x = x(:,1);
t = t0+(0:numel(x)-1)/fs;