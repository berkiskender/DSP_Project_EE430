function spctrgrm(x,fs,t0)
if nargin < 3
    t0 = 0;
end

N = 2.^nextpow2(round(20e-3*fs));   %Window size in samples
win = hamming(N);          %Window
p   = round(N/2);          %Overlap amount in samples
nfft = 2.^nextpow2(N);   %DFT size

[s,f,t] = spectrogram(x,win,p,nfft,fs,'yaxis');
s = 20*log10(abs(s));
t = t + t0;

map = 'default';
clims = [min(min(s(~isinf(s)))) max(max(s))];
imagesc(t,f,s,clims);
set(gca,'YDir','normal');   colormap(map);   colorbar;
title('Spectrogram in dB scale');   xlabel('Time (s)');   ylabel('Frequency (Hz)');
