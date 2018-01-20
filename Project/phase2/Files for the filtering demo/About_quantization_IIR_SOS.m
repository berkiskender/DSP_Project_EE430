% Quantization effect when the second order section (SOS) representation is used.

d = designfilt('lowpassiir', ...
    'FilterOrder', 5, ...
    'PassbandFrequency', 0.1, ...
    'PassbandRipple', 0.1, ...
    'StopbandAttenuation', 100);
s = d.Coefficients();

nbits = 16;
% nbits = 10;
smax = max(max(10*abs(s)));
sq = quantization(s,nbits,smax);
% figure; plot(s(:)); hold on; plot(sq(:),'--');
% fvtool(s,sq);

delta = [1;zeros(1023,1)];   n = 0:numel(delta)-1;
h =  sosfilt(s,delta);
hq = sosfilt(sq,delta);
figure;
subplot(2,1,1); stem(n,h,'filled');  xlabel('n'); title('Impulse response');
subplot(2,1,2); stem(n,hq,'filled'); xlabel('n'); title('Impulse response (with quantized coefficients)');

