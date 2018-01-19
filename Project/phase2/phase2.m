%phase 2-3
t=0:0.00002083:5;
input_sig=100*sin(2*pi*2.*t);



%% 8-bit quantization

encoded_signal_8=uencode(input_sig,8,max(input_sig));

[a,i]=max(input_sig);

decoded_signal_8=udecode(encoded_signal_8,8,max(input_sig));
sample_axis=1:length(input_sig);
%%plot(u,recorded_double);
figure
ax1 = subplot(2,1,1);
stem(ax1,sample_axis,decoded_signal_8);
hold on
stem(ax1,sample_axis,input_sig);
title(ax1,'Original Signal and Decoded Signal for 8 bit quantization');
legend(ax1,'Original Signal','Decoded signal');
error_signal_8=input_sig-decoded_signal_8;
error_signal_power_8=error_signal_8.*error_signal_8;
ax2=subplot(2,1,2);
stem(ax2,sample_axis,error_signal_8);
title(ax2,'Error Signal');
min_error_signal_8=min(error_signal_8);
max_error_signal_8=max(error_signal_8);
mean_error_signal_8=mean(error_signal_power_8);

%% 6-bit quantization

encoded_signal_6=uencode(input_sig,6,max(input_sig));

[a,i]=max(input_sig);

decoded_signal_6=udecode(encoded_signal_6,6,max(input_sig));
sample_axis=1:length(input_sig);
%%plot(u,recorded_double);
figure
ax1 = subplot(2,1,1);
stem(ax1,sample_axis,decoded_signal_6);
hold on
stem(ax1,sample_axis,input_sig);
title(ax1,'Original Signal and Decoded Signal for 6 bit quantization');
legend(ax1,'Original Signal','Decoded Signal');
error_signal_6=input_sig-decoded_signal_6;
error_signal_power_6=error_signal_6.*error_signal_6;
ax2=subplot(2,1,2);
stem(ax2,sample_axis,error_signal_6);
title(ax2,'Error Signal');
min_error_signal_6=min(error_signal_6);
max_error_signal_6=max(error_signal_6);
mean_error_signal_6=mean(error_signal_power_6);
%% 4-bit quantization

encoded_signal_4=uencode(input_sig,4,max(input_sig));

[a,i]=max(input_sig);

decoded_signal_4=udecode(encoded_signal_4,4,max(input_sig));
sample_axis=1:length(input_sig);
%%plot(u,recorded_double);
figure
ax1 = subplot(2,1,1);
stem(ax1,sample_axis,decoded_signal_4);
hold on
stem(ax1,sample_axis,input_sig);
title(ax1,'Original Signal and Decoded Signal for 4 bit quantization');
legend(ax1,'Original Signal','Decoded Signal');
error_signal_4=input_sig-decoded_signal_4;
error_signal_power_4=error_signal_4.*error_signal_4;
ax2=subplot(2,1,2);
stem(ax2,sample_axis,error_signal_4);
title(ax2,'Error Signal');
min_error_signal_4=min(error_signal_4);
max_error_signal_4=max(error_signal_4);
mean_error_signal_4=mean(error_signal_power_4);
