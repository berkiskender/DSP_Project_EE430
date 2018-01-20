% Window-based lowpass FIR filter design.

% Oktay Sipahigil
% 30/12/2016

N = 301;        %Filter length
omgc = 0.37;   %Cut-off frequency (times pi!)
win = rectwin(N)';
scale = @()false;   %If set to true, h[n] is scaled so that H(e^j0)=1.

% Design by using primitive commands
n = 0:N-1;
h = omgc*win.*sinc(omgc*(n-(N-1)/2));
if scale()
    h = h / sum(h);
end

% Design by using fir1
if scale()
    hs = fir1(N-1,omgc,win);
else
    hs = fir1(N-1,omgc,win,'noscale');
end

% Design by using designfilt
d = designfilt('lowpassfir', ...
        'FilterOrder',     N-1, ...
        'CutoffFrequency', omgc, ...
        'DesignMethod',    'window', ...
        'Window',          win, ...
        'ScalePassband',   scale());
hd = d.Coefficients;

% Verify that the alternatives are equal
fprintf('%g\n',sum((hs-h).^2));   assert(sum((hs-h).^2)<eps);
fprintf('%g\n',sum((hd-h).^2));   assert(sum((hd-h).^2)<eps);

% Plot the impulse responses
fvtool(d);
figure; stem(n,h,'filled'); title('h[n]');
xlim([min(n) max(n)]);
xlabel('n');
