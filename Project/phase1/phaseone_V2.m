clear all
close all

%% Sound input from microphone
prompt = 'What is the sampling rate of recorded signal? ';
fs_generated = input(prompt) % Sampling rate of input signal Hz
prompt = 'What is the duration of recorded signal? ';
T_recorded=input(prompt); %input duration in secs
N_recorded=fs_generated*T_recorded; %total number of samples
recorded_sound=audiorecorder(fs_generated,16,1,1); %create an audiorecorder object named recObj

disp('Start speaking.')
recordblocking(recorded_sound, T_recorded);
disp('End of Recording.');

recorded_double = getaudiodata(recorded_sound);
recorded_double=recorded_double';
input_sig=recorded_double;

%% Sound input from file
prompt = 'What is the name of the input audio file? ';
str=input(prompt);
[input_sound,fs_input] = audioread(str); % Take input .mp3 or .wav file as input
input_sound=input_sound(:,1); % Use the single channel of y i1f it has multichannels
input_sig=input_sound;

%% Data generation
prompt = 'Select the type of the signal that you want to generate.\n 0) Sinusoidal signal \n 1) Windowed sinusoidal \n 2) Rectangle Windowed Chirp \n 3) Signal involving multiple components  ';
selection=input(prompt); % Parameter to select which type of signal will be generated

prompt = 'What is the sampling rate of generated signal? ';
fs_generated=input(prompt);% sampling rate of generated signal

%T_generated=; % Total duration of generated signal
%N=fs_generated*T_generated; % total number of samples



if(selection==0)
    
    prompt = 'What is the length of generated sinusoidal? ';
    T_sinusoidal=input(prompt);% length of sine
    
    N_sinusoidal=T_sinusoidal*fs_generated;
    
    %Sinusoidal Signal generation
    prompt = 'What is the amplitude of generated sinusoidal? ';
    A_sinusoidal=input(prompt); % Amplitude of sinusoidal signal

    prompt = 'What is the frequency of generated sinusoidal? ';
    f_sinusoidal=input(prompt); % Frequency of sinusoidal signal
    
    prompt = 'What is the phase of generated sinusoidal? ';
    phase_sinusoidal=input(prompt); % phase of the windowed sinusoidal

    t=linspace(0,T_sinusoidal,N_sinusoidal); % time array involving samples at eact 1/fs secs
    sinusoidal_sig=A_sinusoidal*cos(2*pi*f_sinusoidal.*t+phase_sinusoidal);
    
    input_sig=sinusoidal_sig;
end

    %Windowed sinusoidal signal
if(selection==1)
   
    prompt = 'What is the length of generated windowed sinusoidal? ';
    T_windowed_sinusoidal=input(prompt);% length of windowed sine
    
    N_windowed_sinusoidal=T_windowed_sinusoidal*fs_generated;
    
    prompt = 'What is the starting time of windowed sinusoidal signal? ';
    t0=input(prompt); % Total length was set before, starting time is being set
    
    prompt = 'What is the magnitude of windowed sinusoidal signal? ';
    A_windowed_sinusoidal=input(prompt); % Amplitude of windowed sinusoidal
    
    prompt = 'What is the frequency of windowed sinusoidal signal? ';
    f_windowed_sinusoidal=input(prompt); % Frequency of windowed sinusoidal
    
    prompt = 'What is the phase of windowed sinusoidal signal? ';
    phase_windowed_sinusoidal=input(prompt); % phase of the windowed sinusoidal
    
    t=linspace(t0, T_windowed_sinusoidal+t0, N_windowed_sinusoidal);
    
    windowed_sinusoidal=A_windowed_sinusoidal*cos(2*pi*f_windowed_sinusoidal.*t+phase_windowed_sinusoidal); %rectangular windowed sine is generated, (t-t0 ot t?)

    windowed_sinusoidal=[zeros(1,t0*fs_generated) windowed_sinusoidal];
    
    input_sig=windowed_sinusoidal;

end   

if(selection==2)
    prompt = 'What is the length of generated linear chirp sinusoidal? ';
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
    
    linear_chirp=A_linear_chirp*cos(2*pi*(f_linear_chirp.*t+(t.^2).*bandwidth/(2*T_linear_chirp))+phase_linear_chirp);
    
    linear_chirp=[zeros(1,t0*fs_generated) linear_chirp];
    
    input_sig=linear_chirp;
end
  
% Signal involving multiple components

if(selection==3)
    
    prompt = 'What is the total number of signals? ';
    number_of_signals=input(prompt);
    
    prompt = 'What is the maximum length of generated multiple signals? ';
    T_max=input(prompt);% length total
    total_signal=zeros(1,T_max*fs_generated);
    N_max=T_max*fs_generated;

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

            t=linspace(0,T_max,N_sinusoidal); % time array involving samples at eact 1/fs secs
            sinusoidal_sig=A_sinusoidal*cos(2*pi*f_sinusoidal.*t+phase_sinusoidal);
            
            sinusoidal_sig=[sinusoidal_sig zeros(N_max-length(sinusoidal_sig))];
            
            total_signal=total_signal+sinusoidal_sig;
        end

            %Windowed sinusoidal signal
        if(selection==1)

            prompt = 'What is the length of generated windowed sinusoidal? ';
            T_windowed_sinusoidal=input(prompt);% length of windowed sine

            N_windowed_sinusoidal=T_windowed_sinusoidal*fs_generated;

            prompt = 'What is the starting time of windowed sinusoidal signal? ';
            t0=input(prompt); % Total length was set before, starting time is being set

            prompt = 'What is the magnitude of windowed sinusoidal signal? ';
            A_windowed_sinusoidal=input(prompt); % Amplitude of windowed sinusoidal

            prompt = 'What is the frequency of windowed sinusoidal signal? ';
            f_windowed_sinusoidal=input(prompt); % Frequency of windowed sinusoidal

            prompt = 'What is the phase of windowed sinusoidal signal? ';
            phase_windowed_sinusoidal=input(prompt); % phase of the windowed sinusoidal

            t=linspace(t0, T_windowed_sinusoidal+t0, N_windowed_sinusoidal);

            windowed_sinusoidal=A_windowed_sinusoidal*cos(2*pi*f_windowed_sinusoidal.*t+phase_windowed_sinusoidal); %rectangular windowed sine is generated, (t-t0 ot t?)
            
            windowed_sinusoidal=[zeros(1,t0*fs_generated) windowed_sinusoidal zeros(1,(T_max-t0-T_windowed_sinusoidal)*fs_generated)];
            
            windowed_sinusoidal=[windowed_sinusoidal zeros(N_max-length(windowed_sinusoidal))];
            
            total_signal=total_signal+windowed_sinusoidal;
        end   

        if(selection==2)
            prompt = 'What is the length of generated windowed sinusoidal? ';
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
            
            linear_chirp=A_linear_chirp*cos(2*pi*(f_linear_chirp.*t+(t.^2).*bandwidth/(2*T_linear_chirp))+phase_linear_chirp);
            
            linear_chirp=[zeros(1,t0*fs_generated) linear_chirp zeros(1,(T_max-t0-T_linear_chirp)*fs_generated)];
            
            linear_chirp=[linear_chirp, zeros(N_max-length(linear_chirp))];
            
            total_signal=total_signal+linear_chirp;
        end
        
        
    end
    
    input_sig=total_signal;
end

%% Spectrogram 

prompt = '\nWindow type of STFT:\n 0)Rectangular 1)Hamming 2)Hann 3)Tukey \n4)Triangular 5)Gaussian 6)Blackman 7)Kaiser\n';
type=input(prompt);% window type STFT

prompt = 'Window length of STFT: ';
spec_window_length=input(prompt);% window length of STFT

N_spec=spec_window_length*fs_generated; % total samples in one window

prompt = 'Number of overlapped samples in STFT: ';
overlap=input(prompt);% number of overlapped samples

shift_amount=spec_window_length*fs_generated-overlap; %total amount of shift in samples
l=length(input_sig); %total signal length in samples
total_window=floor(l/(shift_amount)); %total number of possible overlapped or non-overlapping windows in the signal

spectrogram9=zeros(N_spec,total_window); % initiate spectrogram matrix



if(type==0)
    window_type=ones(N_spec,1);
end
    
if(type==1)
    window_type=hamming(N_spec);
end
   
if(type==2)
    window_type=hann(N_spec);
end

if(type==3)
     window_type=tukeywin(N_spec);
end

if(type==4)
     window_type=triang(N_spec);
end
    
if(type==5)
     window_type=gausswin(N_spec);
end

if(type==6)
     window_type=blackman(N_spec);
end

if(type==7)
     window_type=kaiser(N_spec);
end
    
    a=1;
    while (N_spec+(N_spec-overlap)*(total_window-a) > length(input_sig))
        a=a+1;
    end
    
    
    for i=1:(total_window-a);
        signalfft=window_type'.*input_sig((1+(N_spec-overlap)*(i-1)):((N_spec+(N_spec-overlap)*(i-1)))); % one windowed instance of signal to be fed into fft
        spectrogram9(:,i)=(abs(fft(signalfft)));
    end

time=linspace(0,spec_window_length*total_window,total_window); %time matrix for graph
freq=linspace(0,N_spec,N_spec); %frequency matrix for graph
freq=freq'; 

%spectrogram9=spectrogram9'; % transpose of spectrogram

figure,
imagesc(freq,time,abs(spectrogram9));
xlabel('time(s)')
ylabel('freq(Hz)')


%% Benim kodum
   
    signal_kesit=input_sig((1+N_spec*(1-1)):(N_spec*1));
    Fs=fs_generated;
    L=N_spec;
    Y = fft(signal_kesit);
    P2 = abs(Y);
    P2=abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,P1) 
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    [f_row,f_col]=size(P1);
    Den=zeros(f_col,total_window);
    
    Den(:,1)=P1;
    
    a=1;
    while (N_spec+(N_spec-overlap)*(total_window-a) > length(input_sig))
        a=a+1;
    end
    
    for i=2:(total_window-a);
        signal_kesit=window_type'.*input_sig((1+(N_spec-overlap)*(i-1)):((N_spec+(N_spec-overlap)*(i-1))));
        Y = fft(signal_kesit);
        P2 = abs(Y);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        Den(:,i)=P1;
    end
    
 min_Den=min(min(Den));
 max_Den=max(max(Den));
 norm_factor=45/max_Den; % 45 nereden geldi?
 Den=norm_factor*Den;
 y=[0 Fs/2];
 x=[0 T_recorded];

 imagesc(x,y,Den)
 