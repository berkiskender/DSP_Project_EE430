% Quantization effect for the IIR filter

d = designfilt('lowpassiir', ...
    'FilterOrder', 5, ...
    'PassbandFrequency', 0.1, ...
    'PassbandRipple', 0.1, ...
    'StopbandAttenuation', 100);
[b,a] = d.tf();

nbits = 16;
% nbits = 10;
bmax = max(10*abs(b));
amax = max(10*abs(a));
bq = quantization(b,nbits,bmax);
aq = quantization(a,nbits,amax);
% figure; plot(b); hold on; plot(bq,'--');
% figure; plot(a); hold on; plot(aq,'--');

% fvtool(b,a,bq,aq);

delta = [1;zeros(1023,1)];   n = 0:numel(delta)-1;
h =  filter(b,a,delta);
hq = filter(bq,aq,delta);
figure;
subplot(2,1,1); stem(n,h,'filled');  xlabel('n'); title('Impulse response');
subplot(2,1,2); stem(n,hq,'filled'); xlabel('n'); title('Impulse response (with quantized coefficients)');

