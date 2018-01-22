%% Bit Quantization

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
mean_error_signal_8=mean(abs(error_signal_8));

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
mean_error_signal_6=mean(abs(error_signal_6));
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
mean_error_signal_4=mean(abs(error_signal_4));
