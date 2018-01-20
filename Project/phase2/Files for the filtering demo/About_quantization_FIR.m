% Quantization effect on the FIR filters

% x = -110:0.1:110;
% xq = quantization(x,6,100);
% figure; plot(x,x,'--'); hold on; plot(x,xq);
% xlabel('x');   ylabel('x (quantized)');

d = designfilt('lowpassfir', ...
    'FilterOrder', 100, ...
    'PassbandFrequency', 0.1, ...
    'StopbandFrequency', 0.11);
h = d.tf();

nbits = 16;
hmax = max(10*abs(h));
hq = quantization(h,nbits,hmax);
% figure; plot(h); hold on; plot(hq,'--');
% fvtool(h,1,hq,1);

figure;
n = 0:numel(h)-1;
subplot(2,1,1); stem(n,h,'filled');  xlabel('n'); title('Impulse response');
subplot(2,1,2); stem(n,hq,'filled'); xlabel('n'); title('Impulse response (with quantized coefficients)');
