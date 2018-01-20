% Sensor directivity gain pattern measurements.
% Smoothing and interpolation
% The data is circular (angular data).

% Oktay Sipahigil
% 02/01/2017

%% The measurements
x = getdata();   N = numel(x);
tht = 2*(-N/2:N/2-1)/N;   %times pi!
figure; stem(tht,x,'filled'); xlabel('Azimuth angle (\times\pi rad)');
title('Sensor directivity pattern measurements');   ylabel('Gain');

%% Smoothing the measurements
d = designfilt('lowpassfir', ...
    'DesignMethod',    'window', ...
    'FilterOrder',     8, ...
    'CutoffFrequency', 0.4);
fvtool(d);

xx = repmat(x,1,3);
yy = d.filter(xx);
g = mean(d.grpdelay());
y = yy(N+1+g:2*N+g);

figure; stem(tht,x,'filled'); xlabel('Azimuth angle (\times\pi rad)');
title('Sensor directivity pattern measurements');   ylabel('Gain');
hold on; stem(tht,y);
legend('Measurements', 'Smoothed measurements');

%% Obtaining the gain directivity for the in between angles
L = 10;
yy = upsample(y,L);
tt = 2*(-L*N/2:L*N/2-1)/(L*N);
d = designfilt('lowpassfir', ...
    'DesignMethod',    'equiripple', ...
    'FilterOrder', 16*L, ...
    'PassbandFrequency', 0.5/L, ...
    'StopbandFrequency', 1.5/L);
fvtool(d);
yy = repmat(yy,1,3);
z = L*d.filter(yy);
g = mean(d.grpdelay());
z = z(L*N+1+g:L*2*N+g);
figure; stem(tht,x,'filled'); xlabel('Azimuth angle (\times\pi rad)');
title('Sensor directivity pattern measurements');   ylabel('Gain');
hold on; plot(tt,z,'Linewidth',2);
legend('Measurements', 'Smoothed and then interpolated');
