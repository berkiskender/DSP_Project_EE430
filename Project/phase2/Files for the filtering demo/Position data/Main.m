% A noisy data and its approximate derivative

% Oktay Sipahigil
% 02/01/2017

%% The noisy position data
[x,fs] = getdata();
t = (0:numel(x)-1)/fs;
figure;
plot(t,x); xlabel('Time (s)'); ylabel('Position (m)');

%% Obtain the speed (1)
z = [0 diff(x)];
figure;
subplot(2,1,1); plot(t,x); xlabel('Time (s)'); ylabel('Position (m)');
subplot(2,1,2); plot(t,z); xlabel('Time (s)'); ylabel('Speed (m/s)');
xlim([0 3]);

%% Obtain the speed (2)
d = designfilt('differentiatorfir', ...
    'DesignMethod', 'equiripple', ...
    'FilterOrder', 50, ...
    'PassbandFrequency', 5, ...
    'StopbandFrequency', 7, ...
    'SampleRate', fs);
fvtool(d);

y = d.filter(x) * fs;
g = mean(d.grpdelay());
figure;
subplot(2,1,1); plot(t,x); xlabel('Time (s)'); ylabel('Position (m)');
subplot(2,1,2); plot(t(1+g:end),y(1+g:end)); xlabel('Time (s)'); ylabel('Speed (m/s)');
xlim([0 3]);

