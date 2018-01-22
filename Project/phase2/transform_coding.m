
%% Transform Coding

prompt='Do you want to analyze a \n0)Recorded sound \n1)Sound data from a file \n2)Generated data ?';
choice = input(prompt);
%% Sound input from microphone
if (choice==0)
    prompt = 'What is the sampling rate of recorded signal? ';
    fs_generated = input(prompt); % Sampling rate of input signal Hz
    prompt = 'What is the duration of recorded signal? ';
    T_signal=input(prompt); %input duration in secs
    N=fs_generated*T_signal; %total number of samples
    recorded_sound=audiorecorder(fs_generated,16,1,1); %create an audiorecorder object named recObj

    disp('Start speaking.')
    recordblocking(recorded_sound, T_signal);
    disp('End of Recording.');

    recorded_double = getaudiodata(recorded_sound);
    recorded_double=recorded_double';
    input_sig=recorded_double;
end
%% Sound input from file
if(choice==1)
    prompt = 'What is the name of the input audio file? ';
    str=input(prompt);
    [input_sound,fs_input] = audioread(str); % Take input .mp3 or .wav file as input
    input_sound=input_sound(:,1); % Use the single channel of y i1f it has multichannels
    input_sig=input_sound';
    fs_generated=fs_input;
    T_signal=length(input_sig)/fs_generated;
    N=length(input_sig);
end
%% Data generation
if(choice==2)
    prompt = 'Select the type of the signal that you want to generate.\n 0) Sinusoidal signal \n 1) Windowed sinusoidal \n 2) Rectangle Windowed Chirp \n 3) Signal involving multiple components  ';
    selection=input(prompt); % Parameter to select which type of signal will be generated

    prompt = 'What is the sampling rate of generated signal? ';
    fs_generated=input(prompt);% sampling rate of generated signal

    %T_generated=; % Total duration of generated signal
    %N=fs_generated*T_generated; % total number of samples



    if(selection==0)

        prompt = 'What is the length of generated sinusoidal? ';
        T_signal=input(prompt);% length of sine


        N_signal=T_signal*fs_generated;

        %Sinusoidal Signal generation
        prompt = 'What is the amplitude of generated sinusoidal? ';
        A_sinusoidal=input(prompt); % Amplitude of sinusoidal signal

        prompt = 'What is the frequency of generated sinusoidal? ';
        f_sinusoidal=input(prompt); % Frequency of sinusoidal signal

        prompt = 'What is the phase of generated sinusoidal? ';
        phase_sinusoidal=input(prompt); % phase of the windowed sinusoidal

        t=linspace(0,T_signal,N_signal); % time array involving samples at eact 1/fs secs
        sinusoidal_sig=A_sinusoidal*cos(2*pi*f_sinusoidal.*t+phase_sinusoidal);

        input_sig=sinusoidal_sig;
    end

        %Windowed sinusoidal signal
    if(selection==1)

        prompt = 'What is the length of generated windowed sinusoidal? ';
        T_signal=input(prompt);% length of windowed sine

        N_signal=T_signal*fs_generated;

        prompt = 'What is the starting time of windowed sinusoidal signal? ';
        t0=input(prompt); % Total length was set before, starting time is being set

        prompt = 'What is the magnitude of windowed sinusoidal signal? ';
        A_windowed_sinusoidal=input(prompt); % Amplitude of windowed sinusoidal

        prompt = 'What is the frequency of windowed sinusoidal signal? ';
        f_windowed_sinusoidal=input(prompt); % Frequency of windowed sinusoidal

        prompt = 'What is the phase of windowed sinusoidal signal? ';
        phase_windowed_sinusoidal=input(prompt); % phase of the windowed sinusoidal

        t=linspace(t0, T_signal+t0, N_signal);

        windowed_sinusoidal=A_windowed_sinusoidal*cos(2*pi*f_windowed_sinusoidal.*t+phase_windowed_sinusoidal); %rectangular windowed sine is generated, (t-t0 ot t?)

        windowed_sinusoidal=[zeros(1,t0*fs_generated) windowed_sinusoidal];

        input_sig=windowed_sinusoidal;

        T_signal=t0+T_signal;


    end   

    if(selection==2)
        prompt = 'What is the length of generated linear chirp sinusoidal? ';
        T_signal=input(prompt);% length of chirp



        prompt = 'What is the starting time of linear chirp? ';
        t0=input(prompt);% starting time of chirp

        N_signal=(T_signal-t0)*fs_generated;

        prompt = 'What is the amplitude of linear chirp? ';
        A_linear_chirp=input(prompt);% amplitude of the chirp

        prompt = 'What is the frequency of linear chirp? ';
        f_linear_chirp=input(prompt);% frequency of the chirp

        prompt = 'What is the phase of linear chirp? ';
        phase_linear_chirp=input(prompt);% phase of the chirp

        prompt = 'What is the bandwidth of linear chirp? ';
        bandwidth=input(prompt);% bandwidth of the chirp

        t=linspace(t0, T_signal, N_signal);

        linear_chirp=A_linear_chirp*cos(2*pi*(f_linear_chirp.*t+(t.^2).*bandwidth/(2*(T_signal-t0)))+phase_linear_chirp);

        linear_chirp=[zeros(1,t0*fs_generated) linear_chirp];

        input_sig=linear_chirp;
    end

    % Signal involving multiple components

    if(selection==3)

        prompt = 'What is the total number of signals? ';
        number_of_signals=input(prompt);

        prompt = 'What is the maximum length of generated multiple signals? ';
        T_signal=input(prompt);% length total
        input_sig=zeros(1,T_signal*fs_generated);
        N_signal=T_signal*fs_generated;

        for i=1:number_of_signals


            fprintf('\n\n %d. signal\n',i)
            prompt = 'Select the type of the signal that you want to generate.\n 0) Sinusoidal signal \n 1) Windowed sinusoidal \n 2) Rectangle Windowed Chirp \n  ';
            selection=input(prompt); % Parameter to select which type of signal will be generated

            if(selection==0)

                prompt = 'What is the length of generated sinusoidal? ';
                T_max_sinusoidal=input(prompt);% length of sine

                N_sinusoidal=T_max_sinusoidal*fs_generated;

                %Sinusoidal Signal generation
                prompt = 'What is the amplitude of generated sinusoidal? ';
                A_sinusoidal=input(prompt); % Amplitude of sinusoidal signal

                prompt = 'What is the frequency of generated sinusoidal? ';
                f_sinusoidal=input(prompt); % Frequency of sinusoidal signal

                prompt = 'What is the phase of generated sinusoidal? ';
                phase_sinusoidal=input(prompt); % phase of the windowed sinusoidal

                t=linspace(0,T_max_sinusoidal,N_sinusoidal); % time array involving samples at eact 1/fs secs
                sinusoidal_sig=A_sinusoidal*cos(2*pi*f_sinusoidal.*t+phase_sinusoidal);

                sinusoidal_sig=[sinusoidal_sig zeros(1,N_signal-length(sinusoidal_sig))];

                input_sig=input_sig+sinusoidal_sig;


            end

                %Windowed sinusoidal signal
            if(selection==1)

                prompt = 'What is the length of generated windowed sinusoidal? ';
                T_windowed_sinusoidal=input(prompt);% length of windowed sine

                N_signal=T_windowed_sinusoidal*fs_generated;

                prompt = 'What is the starting time of windowed sinusoidal signal? ';
                t0=input(prompt); % Total length was set before, starting time is being set

                prompt = 'What is the magnitude of windowed sinusoidal signal? ';
                A_windowed_sinusoidal=input(prompt); % Amplitude of windowed sinusoidal

                prompt = 'What is the frequency of windowed sinusoidal signal? ';
                f_windowed_sinusoidal=input(prompt); % Frequency of windowed sinusoidal

                prompt = 'What is the phase of windowed sinusoidal signal? ';
                phase_windowed_sinusoidal=input(prompt); % phase of the windowed sinusoidal

                t=linspace(t0, T_windowed_sinusoidal+t0, N_signal);

                windowed_sinusoidal=A_windowed_sinusoidal*cos(2*pi*f_windowed_sinusoidal.*t+phase_windowed_sinusoidal); %rectangular windowed sine is generated, (t-t0 ot t?)

                windowed_sinusoidal=[zeros(1,t0*fs_generated) windowed_sinusoidal zeros(1,(T_signal-t0-T_windowed_sinusoidal)*fs_generated)];

                windowed_sinusoidal=[windowed_sinusoidal zeros(1,N_signal-length(windowed_sinusoidal))];

                input_sig=input_sig+windowed_sinusoidal;

                T_signal=t0+T_signal;
            end   

            if(selection==2)
                prompt = 'What is the length of generated windowed linear chirp? ';
                T_linear_chirp=input(prompt);% length of chirp

                N_linear_chirp=T_linear_chirp*fs_generated;

                prompt = 'What is the starting time of linear chirp? ';
                t0=input(prompt);% starting time of chirp

                prompt = 'What is the amplitude of linear chirp? ';
                A_linear_chirp=input(prompt);% amplitude of the chirp

                prompt = 'What is the frequency of linear chirp? ';
                f_linear_chirp=input(prompt);% frequency of the chirp

                prompt = 'What is the phase of linear chirp? ';
                phase_linear_chirp=input(prompt);% phase of the chirp

                prompt = 'What is the bandwidth of linear chirp? ';
                bandwidth=input(prompt);% bandwidth of the chirp

                t=linspace(t0, T_linear_chirp+t0, N_linear_chirp);

                linear_chirp=A_linear_chirp*cos(2*pi*(f_linear_chirp.*(t-t0)+((t-t0).^2).*bandwidth/(2*T_linear_chirp))+phase_linear_chirp);

                linear_chirp=[zeros(1,t0*fs_generated) linear_chirp zeros(1,(T_signal-t0-T_linear_chirp)*fs_generated)];

                linear_chirp=[linear_chirp, zeros(1,N_signal-length(linear_chirp))];

                input_sig=input_sig+linear_chirp;
            end


        end

        input_sig=input_sig;
    end
N=N_signal;    
end

%% Input signal spectrogram

spectrogram_group9(input_sig, fs_generated, length(input_sig)/fs_generated);

%% DFT
%Coding (Supply desired number of nonzero coefficients)
%Sort absolute values of the coefficients and then determine the
%corresponding threshold automatically for each input audio signal

%Listen, comment on quality, plot spectrogram of input and output, compute error signal, average error 
partition_amount_dft=500;

% Sorting coeffs with respect to absolute values
% sorted_input_sig_fft=sort(abs(fft(input_sig)));
% Selecting threshold value (Decides on compression amount) (ascending order sorted so value represents the amount of eliminated elements)
% threshold=sorted_input_sig_fft(ceil(compression_amount));
% max_sorted_input=max(sorted_input_sig_fft);

% Using DFT (Partition the input signal and apply quantization to frequency domain coefficients)
N_blocksize=floor(length(input_sig)/partition_amount_dft); % block size to be used in DFT
compression_amount=9000*N_blocksize/10000;

input_sig_buffered=buffer(input_sig ,N_blocksize); % Partition the signal

input_sig_buffered_fft=fft(input_sig_buffered); % take fft of each buffered partition

for i=1:partition_amount_dft
    thresholded_partition=input_sig_buffered_fft(:,i);
    sorted_input_sig_fft=sort(abs(thresholded_partition));
    threshold=sorted_input_sig_fft(ceil(compression_amount));
    thresholded_partition(abs(thresholded_partition) < threshold)=0;
    input_sig_buffered_fft(:,i)=thresholded_partition;
end
% input_sig_buffered_fft(abs(input_sig_buffered_fft) < threshold)=0; % make values zero if they have abs value smaller than threshold

input_sig_buffered_thresholded=ifft(input_sig_buffered_fft);
input_sig_thresholded=input_sig_buffered_thresholded(:); % obtain a vector from buffered & thresholded matrix

% Error signal computation
input_sig_thresholded=input_sig_thresholded(1:length(input_sig));
err_dft_comp=input_sig-input_sig_thresholded';

% Average power computation
err_avg_pwr=err_dft_comp.^2;
err_avg_pwr=mean(err_avg_pwr);
err_tot_pwr=sum(err_dft_comp.^2);
input_sig_avg_pwr=mean(input_sig.^2);


input_sig_thresholded=input_sig_thresholded';
avg_pwr_input_sig_thresholded=mean(input_sig_thresholded.^2);

numzeros_dft=nnz(~real(input_sig_buffered_fft));

%% Figure plots FFT

% err signal
figure, 
plot(err_dft_comp);
ylabel('magnitude')
xlabel('samples')
str=sprintf('DFT error signal, average power of error: %f, \n average power of signal %f, \n partition amount: %d', err_avg_pwr,input_sig_avg_pwr, partition_amount_dft);
title(str);

% FFT of input sig
figure, 
subplot(1,2,1);
plot(fftshift(abs(fft(input_sig))));
title('original input signal DFT')
ylabel('magnitude')
xlabel('samples')
% FFT of thresholded sig
subplot(1,2,2); 
plot(abs(fftshift(fft(input_sig_thresholded))));
ylabel('magnitude')
xlabel('samples')
str=sprintf('DFT thresholded, percentage of deleted coefficients: %f, \n partition amount: %d',100*(compression_amount/N_blocksize), partition_amount_dft);
title(str);

%% spectrogram of compressed signal
fs_generated=44100;
spectrogram_group9(input_sig_thresholded, fs_generated, length(input_sig)/fs_generated);

%% Play compressed and original using DFT

% input_sig_thresholded=ifft(input_sig_thresholded_fft);
% sound(real(input_sig_thresholded_dct), fs_generated);
% sound(real(input_sig_thresholded), fs_generated);
% sound(real(input_sig), fs_generated);

%% DCT

%Listen, comment on quality, plot spectrogram of input and output, compute error signal, average error 
partition_amount_dct=500;

% Using DCT (Partition the input signal and apply quantization to frequency domain coefficients)
N_blocksize_dct=floor(length(input_sig)/partition_amount_dct); % block size to be used in DCT
compression_amount=9000*N_blocksize_dct/10000;
% dct
input_sig_buffered=buffer(input_sig ,N_blocksize_dct); % Partition the signal
% audioread
input_sig_buffered_dct=dct(input_sig_buffered,[],1,'Type',2); % take dct of each buffered partition

for i=1:partition_amount_dct
    thresholded_partition=input_sig_buffered_dct(:,i);
    sorted_input_sig_dct=sort(abs(thresholded_partition));
    threshold_dct=sorted_input_sig_dct(ceil(compression_amount));
    thresholded_partition(abs(thresholded_partition) < threshold_dct)=0;
    input_sig_buffered_dct(:,i)=thresholded_partition;
end
% input_sig_buffered_fft(abs(input_sig_buffered_fft) < threshold)=0; % make values zero if they have abs value smaller than threshold

input_sig_buffered_thresholded=idct(input_sig_buffered_dct,[],1,'Type',2);
input_sig_thresholded=input_sig_buffered_thresholded(:); % obtain a vector from buffered & thresholded matrix

% Error signal computation
input_sig_thresholded_dct=input_sig_thresholded(1:length(input_sig));
err_dct_comp=input_sig-input_sig_thresholded_dct';

% Average power computation
err_avg_pwr=err_dct_comp.^2;
err_avg_pwr=mean(err_avg_pwr);
err_tot_pwr=sum(err_dct_comp.^2);
input_sig_avg_pwr=mean(input_sig.^2);

input_sig_thresholded_dct=input_sig_thresholded_dct';
avg_pwr_input_sig_thresholded_dct=mean(input_sig_thresholded_dct.^2);
numzeros_dct=nnz(~real(input_sig_buffered_dct));
%% Figure plots DCT

% err signal
figure, 
plot(err_dct_comp);
ylabel('magnitude')
xlabel('samples')
str=sprintf('DCT-2 error signal, average power of error: %f, \n average power of signal %f, \n partition amount: %d', err_avg_pwr,input_sig_avg_pwr, partition_amount_dct);
title(str);

% FFT of input sig
figure, 
subplot(1,2,1);
plot(abs(fftshift(fft(input_sig))));
title('original input signal DCT')
ylabel('magnitude')
xlabel('samples')
% FFT of thresholded sig
subplot(1,2,2); 
plot(abs(fftshift(fft(input_sig_thresholded))));
ylabel('magnitude')
xlabel('samples')
str=sprintf('DCT-2 thresholded, percentage of deleted coefficients: %f, \n partition amount: %d',100*(compression_amount/N_blocksize_dct),partition_amount_dct);
title(str);

%% MDCT

% partition input signal into overlapping blocks of 50%
%Listen, comment on quality, plot spectrogram of input and output, compute error signal, average error 
partition_amount_mdct=500;

% Using MDCT (Partition the input signal and apply quantization to frequency domain coefficients)
N2_blocksize_mdct=floor(length(input_sig)/(partition_amount_mdct/2)); % block size to be used in DCT

N_blocksize_mdct=round(N2_blocksize_mdct/2);
compression_amount=9000*N2_blocksize_mdct/10000;
%=======
N_blocksize_mdct=floor(N2_blocksize_mdct/2);
compression_amount=9000*N_blocksize_mdct/10000;
%>>>>>>> 0537f7d4b2996cf39904891c7fee1cac37400979
%=======

N_blocksize_mdct=round(N2_blocksize_mdct/2);
compression_amount=9000*N2_blocksize_mdct/10000;

N_blocksize_mdct=floor(N2_blocksize_mdct/2);
compression_amount=9000*N_blocksize_mdct/10000;

%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
% mdct
input_sig_buffered=buffer(input_sig ,N2_blocksize_mdct,round(N2_blocksize_mdct/2)); % Partition the signal

y_input_sig_buffered=zeros(N_blocksize_mdct,partition_amount_mdct+1);

%<<<<<<< HEAD
%<<<<<<< HEAD
%=======
%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
 for k=1:partition_amount_mdct+1
 for i=1:N_blocksize_mdct
     if(i<(floor(N_blocksize_mdct/2)+1))
         y_input_sig_buffered(i,k)=-input_sig_buffered(floor(i+3*N_blocksize_mdct/2),k)-input_sig_buffered(floor(3*N_blocksize_mdct/2-1-i),k);
    else
%<<<<<<< HEAD
        y_input_sig_buffered(i,k)=input_sig_buffered(floor(i-N_blocksize_mdct/2),k)-input_sig_buffered(floor(3*N_blocksize_mdct/2-1-i),k);
     end
 end
 end
%=======
%=======
        y_input_sig_buffered(i,k)=input_sig_buffered(floor(i-N_blocksize_mdct/2),k)-input_sig_buffered(floor(3*N_blocksize_mdct/2-1-i),k); 
%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
for k=1:partition_amount_mdct
for i=0:N_blocksize_mdct-1
    if(i<floor(N_blocksize_mdct/2)+1)
        y_input_sig_buffered(i+1,k)=-input_sig_buffered(floor(i+3*N_blocksize_mdct/2)+1,k)-input_sig_buffered(floor(3*N_blocksize_mdct/2-1-i)+1,k);
    else
        y_input_sig_buffered(i+1,k)=input_sig_buffered(floor(i-N_blocksize_mdct/2)+1,k)-input_sig_buffered(floor(3*N_blocksize_mdct/2-1-i)+1,k); 
%<<<<<<< HEAD
%>>>>>>> 0537f7d4b2996cf39904891c7fee1cac37400979
%=======
%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
    end
    end
end

%<<<<<<< HEAD
%<<<<<<< HEAD
%=======
%=======
%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
input_sig_buffered_mdct=dct(y_input_sig_buffered,[],1,'Type',4);

for i=1:partition_amount_mdct
    thresholded_partition=input_sig_buffered_mdct(:,i);
    sorted_input_sig_mdct=sort(abs(thresholded_partition));
    threshold_mdct=sorted_input_sig_mdct(ceil(compression_amount));
    thresholded_partition(abs(thresholded_partition) < threshold_mdct)=0;
    input_sig_buffered_mdct(:,i)=thresholded_partition;
end

% IDCT
Beta=idct(input_sig_buffered_mdct,[],1,'Type',4); 

alpha=zeros(1,N2_blocksize_mdct);

for n=0:floor(N_blocksize_mdct/2)-1
    alpha(n+1)=Beta((n+1)+floor(N_blocksize_mdct/2));
end

%!!!
for n=(floor(N_blocksize_mdct/2)):(floor(3*N_blocksize_mdct/2)-1)
    alpha(n+1)=(-1)*Beta(ceil(3*N_blocksize_mdct/2)-(n+1));
end

for n=(floor(3*N_blocksize_mdct/2)):(floor(N2_blocksize_mdct)-1)
    alpha(n+1)=-1*Beta(-1*floor(3*N_blocksize_mdct/2)+(n+1));
end

xr=zeros(N_blocksize_mdct);



%% 
for i=0:length(a)/2-1
if(i<floor(length(a)/4))
b(i+1)=-a(floor(i+3*length(a)/4)+1)-a(floor(3*length(a)/4-1-(i))+1);
else
b(i+1)=a(floor(i-length(a)/4)+1)-a(floor(3*length(a)/4-1-(i))+1);
end
end




% input_sig_buffered_mdct=dct(input_sig_buffered,[],1,'Type',4); % take dct of each buffered partition

% for i=1:partition_amount_dct
%     thresholded_partition=input_sig_buffered_mdct(:,i);
%     sorted_input_sig_dct=sort(abs(thresholded_partition));
%     threshold_dct=sorted_input_sig_dct(ceil(compression_amount));
%     thresholded_partition(abs(thresholded_partition) < threshold_dct)=0;
%     input_sig_buffered_mdct(:,i)=thresholded_partition;
% end

%<<<<<<< HEAD



%>>>>>>> 0537f7d4b2996cf39904891c7fee1cac37400979

%=======
%>>>>>>> b75b29523feda09decf7596d49b6cd970bd40dfe
