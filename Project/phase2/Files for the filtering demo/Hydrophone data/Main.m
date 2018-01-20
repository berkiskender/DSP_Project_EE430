% A hydrophone data example
% Applications of high pass filtering, notch filtering.
% Envelope finding is through Hilbert transform.

% Oktay Sipahigil
% 02/01/2017

%% The data
[x,fs] = getdata();
t = (0:numel(x)-1)/fs;

figure;   plot(t,x);
xlabel('Time (sec)');   ylabel('Voltage (mV)');   title('Hydrophone data');
figure; spctrgrm(x,fs);

% nfft = 2^nextpow2(numel(x));
% X = fft(x,nfft);   X = X(1:nfft/2+1);
% f = (0:nfft/2)*fs/nfft;
% figure; semilogx(f,10*log10(X.*conj(X))); grid minor;
% xlabel('Frequency (Hz)');   ylabel('|X(e^{j\omega})| in dB');   title('Magnitude characteristics of hydrophone data');


%% High-pass filtering
dhp = designfilt('highpassiir', ...
    'DesignMethod',        'butter', ...
    'FilterOrder',         6, ...
    'HalfPowerFrequency',  10/fs*2);
fvtool(dhp);
y = dhp.filter(x);
figure;
subplot(2,1,1);   plot(t,x);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Hydrophone data');
subplot(2,1,2);   plot(t,y);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Highpass filtered hydrophone data');
figure; spctrgrm(y,fs);

%% Band-stop filtering
dn = designfilt('bandstopiir', ...
    'DesignMethod',        'butter', ...
    'FilterOrder',         6, ...
    'HalfPowerFrequency1', 0.9e3/fs*2, ...
    'HalfPowerFrequency2', 1e3/fs*2);
fvtool(dn);
z = dn.filter(y);
figure;
subplot(3,1,1);   plot(t,x);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Hydrophone data');
subplot(3,1,2);   plot(t,y);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Highpass filtered hydrophone data');
subplot(3,1,3);   plot(t,z);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Highpass and notch filtered hydrophone data');
figure; spctrgrm(z,fs);

%% Hilbert transform filtering
dh = designfilt('hilbertfir', ...
    'DesignMethod',        'equiripple', ...
    'FilterOrder',         200, ...
    'TransitionWidth',     100/fs*2);
fvtool(dh);
% xx = cos(2*pi*1e3*t);
% yy = dh.filter(xx);
% ty = t - floor(mean(dh.grpdelay()))/fs;
% figure; plot(t,xx); hold on; plot(ty,yy);

w = dh.filter(z);
g = floor(mean(dh.grpdelay()));
n0 = 250;
tt = t(n0+1:end-g);
zz = z(1:end-g);   zz = zz(n0+1:end);
ww = w(g+1:end);   ww = ww(n0+1:end);
figure; plot(zz); hold on; plot(ww);
e = abs(zz + 1j*ww);
figure;
subplot(3,1,1);   plot(t,x);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Hydrophone data');
subplot(3,1,2);   plot(t,y);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Highpass filtered hydrophone data');
subplot(3,1,3);   plot(t,z);   xlabel('Time (s)');   ylabel('Voltage (mV)');   title('Highpass and notch filtered hydrophone data');
hold on;   plot(tt,[1;-1]*e,'Linewidth',2);
