
% h=firpm(511,[0 .003906 .015625 .5]*2,[2 2 0 0]);
h=prototype_filter;
figure, plot(h);
title('h[n]')
xlabel('samples')
ylabel('magnitude')

x=-pi+pi/256:pi/256:pi;

figure, plot(x,abs(fftshift(fft(h))));
title('Magnitude response of h[n]')
xlabel('frequency (rad)')
ylabel('magnitude')

figure, plot(x,phase(fftshift(fft(h))));
title('Phase response of h[n]')
xlabel('frequency (rad)')
ylabel('phase')
length(h);

for k=0:31
    for n=1:512
        h32(n,k+1)=h(n)*cos((pi/64)*(2*k+1)*(n-17));
    end
end

figure, plot(x,abs(fftshift(fft(h32))));
title('Magnitude response of h0 to h31')
xlabel('frequency (rad)')
ylabel('magnitude')

figure, plot(x,phase(fftshift(fft(h32(:,3)))));
title('Phase response of h0 to h31')
xlabel('frequency (rad)')
ylabel('phase')


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

x128=-pi+pi/(0.5*length(singing_44100_128)):pi/(0.5*length(singing_44100_128)):pi;
x96=-pi+pi/(0.5*length(singing_44100_96)):pi/(0.5*length(singing_44100_96)):pi;
x192=-pi+pi/(0.5*length(singing_44100_192)):pi/(0.5*length(singing_44100_192)):pi;
xorig=-pi+pi/(0.5*length(singing_44100_orig)):pi/(0.5*length(singing_44100_orig)):pi;;

x128=x128';
x96=x96';
x192=x192';
xorig=xorig';

figure,
subplot(1,2,1);
plot(xorig,fftshift(abs(fft(singing_44100_orig))));
title('Original Signal')
ylabel('Magnitude')
xlabel('frequency (rad) (min-max=-+22050 Hz)')
subplot(1,2,2);
plot(x96,fftshift(abs(fft(singing_44100_96))));
title('compressed signal rate=96 kbps')
ylabel('Magnitude')
xlabel('frequency (rad) (min-max=-+22050 Hz)')
figure,
subplot(1,2,1);
plot(x128,fftshift(abs(fft(singing_44100_128))));
title('compressed signal rate=128 kbps')
ylabel('Magnitude')
xlabel('frequency (rad) (min-max=-+22050 Hz)')
subplot(1,2,2);
plot(x192,fftshift(abs(fft(singing_44100_192))));
title('compressed signal rate=192 kbps')
ylabel('Magnitude')
xlabel('frequency (rad) (min-max=-+22050 Hz)')

