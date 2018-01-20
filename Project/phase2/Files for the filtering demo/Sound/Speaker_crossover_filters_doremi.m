% Split full-range audio into three channels with low, mid and high frequency bands.
% Each channel may be sent to different speakers (Woofer, midrange driver and tweeter)
% See the section titled "Speaker Crossover Filters" in Signal Processing Toolbox User's Guide.
%
% To stop playing audio:
%    clear sound 

% Oktay Sipahigil
% 02/01/2017

fs = 16000;   %Sampling frequency in Hz
dur = 10;

t = (0:dur*fs-1)/fs;
x = doremi(t) + 0.001*randn(size(t));
figure; spctrgrm(x,fs,t(1));
% sound(x,fs);

%% Filter desgins
omg0 = 1088/fs*2;   %First delimiter frequency (times pi!)
omg1 = 2536/fs*2;   %Second delimiter frequency (times pi!)
ord = 6;   %Filter order
% ord = 20;   %Filter order
rp = 1;   %Passband ripple in dB

dlp = designfilt('lowpassiir', ...
    'DesignMethod',      'cheby1', ...
    'FilterOrder',       ord, ...
    'PassbandFrequency', omg0, ...
    'PassbandRipple',    rp);

% dlp = designfilt('lowpassiir', ...
%     'DesignMethod',      'ellip', ...
%     'FilterOrder',       ord, ...
%     'PassbandFrequency', omg0, ...
%     'PassbandRipple',    rp, ...
%     'StopbandAttenuation', 100);

dbp = designfilt('bandpassiir', ...
    'DesignMethod',       'cheby1', ...
    'FilterOrder',        ord, ...
    'PassbandFrequency1', omg0, ...
    'PassbandFrequency2', omg1, ...
    'PassbandRipple',     rp);

dhp = designfilt('highpassiir', ...
    'DesignMethod',      'cheby1', ...
    'FilterOrder',       ord, ...
    'PassbandFrequency', omg1, ...
    'PassbandRipple',    rp);

fvtool(dlp,dbp,dhp);
% fvtool(dlp);
% fvtool(dbp);
% fvtool(dhp);

%% Filter data
ylp = dlp.filter(x);
ybp = dbp.filter(x);
yhp = dhp.filter(x);
figure; spctrgrm(ylp,fs); title('Lowpass channel  (Woofer)');
figure; spctrgrm(ybp,fs); title('Bandpass channel (Midrange)');
figure; spctrgrm(yhp,fs); title('Highpass channel (Tweeter)');
% sound(ylp,fs);
% sound(ybp,fs);
% sound(yhp,fs);

%% Equivalent design by calling the cheby1 function
% [sb0,sa0] = cheby1(ord,rp,omg0);
% [sb1,sa1] = cheby1(ord/2,rp,[omg0 omg1]);
% [sb2,sa2] = cheby1(ord,rp,omg1,'high');
% [blp,alp] = dlp.tf();   assert(sum((blp-sb0).^2)<eps);   assert(sum((alp-sa0).^2)<eps);
% [bbp,abp] = dbp.tf();   assert(sum((bbp-sb1).^2)<eps);   assert(sum((abp-sa1).^2)<eps);
% [bhp,ahp] = dhp.tf();   assert(sum((bhp-sb2).^2)<eps);   assert(sum((ahp-sa2).^2)<eps);
