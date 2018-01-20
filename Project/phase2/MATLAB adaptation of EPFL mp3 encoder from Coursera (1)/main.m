
% h=firpm(511,[0 .003906 .015625 .5]*2,[2 2 0 0]);
h=prototype_filter;
figure, plot(h);
title('h[n]')
xlabel('samples')
ylabel('magnitude')

figure, plot(abs(fft(h)));
title('Magnitude response of h[n]')
xlabel('frequency')
ylabel('magnitude')

figure, plot(phase(fft(h)));
title('Phase response of h[n]')
xlabel('frequency')
ylabel('phase')
length(h);

for k=0:31
    for n=1:512
        h32(n,k+1)=h(n)*cos((pi/64)*(2*k+1)*(n-17));
    end
end

figure, plot(abs(fft(h32)));
title('Magnitude response of h0 to h31')
xlabel('frequency')
ylabel('magnitude')

t=linspace(0,1,512);
t=t';
x=cos(2*pi*20*t); % example sinusoidal
h=prototype_filter;
s = subband_filtering(x,h);
X = scaled_fft_db(x)
figure, plot(X);
title ('257-point array of spectral magnitude values in dBs, normalized to 96dB');
ylabel('magnitude (dB)')

Encoder_main_script;
