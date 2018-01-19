function [y1,yN] = spectrogram_group9(input_sig, fs_generated, T_signal)
%% Spectrogram 

prompt = '\nWindow type of STFT:\n 0)Rectangular 1)Hamming 2)Hann 3)Tukey \n4)Triangular 5)Gaussian 6)Blackman 7)Kaiser\n';
type=input(prompt);% window type STFT

prompt = 'Window length of STFT: ';
spec_window_length=input(prompt);% window length of STFT

N_spec=floor(spec_window_length*fs_generated); % total samples in one window

prompt = 'Overlap length in STFT: ';
overlap=floor(input(prompt)*fs_generated);% number of overlapped samples

shift_amount=spec_window_length*fs_generated-overlap; %total amount of shift in samples
l=length(input_sig); %total signal length in samples
total_window=floor(l/(shift_amount)); %total number of possible overlapped or non-overlapping windows in the signal

spectrogram9=zeros((N_spec),total_window); % initiate spectrogram matrix



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

%figure,
%imagesc(freq,time,abs(spectrogram9));
%xlabel('time(s)')
%ylabel('freq(Hz)')


%% Plot and calculation of STFT
   
    signal_kesit=input_sig((1+N_spec*(1-1)):(N_spec*1));
    Fs=fs_generated;
    L=N_spec;
    Y = fft(signal_kesit);
    P2 = abs(Y);
    P2=10*log10(abs(Y/L).^2);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
%     plot(f,P1) 
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
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
        %P2=10.*log10(abs(Y/L).^2); %normalization
        P2=10.*log10(abs(Y).^2);  % NOT normalized
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        Den(:,i)=P1;
    end
    
 y=[0 Fs/2];
 x=[0 T_signal];
 
min_den=min(Den(:));
max_den=max(Den(:));

%Den(Den<0.1) = 0;

clims=[min_den max_den] ;
%im=imagesc(x,y,10*log10(Den.^2),clims);
figure,imagesc(x,y,Den,clims);
title('spectrogram of input signal')
xlabel('time(s)');
ylabel('freq(Hz)');
colorbar;

end