function X = scaled_fft_db(x)
% takes one input and returns 251x1 array of spectrum estimate
%Compute Hann window
w1=hann(512);
w2=zeros(512,1);
alpha=sqrt(511/(sum(w1.^2)));
for n=1:512
w2(n)=(alpha/2)*(1-cos((2*pi*n)/(512-1)));
end

% windowing (time domain mult.)
windowed_x=x.*w2;
% Compute fft
windowed_X_fft = fft(windowed_x);
% Normalize by length (512)
windowed_X_fft = windowed_X_fft/length(windowed_X_fft);
% discard after 257th element since real
windowed_X_fft =  windowed_X_fft(1:257);
% Convert magnitude to dB, if original is zero, set to -100 dB
windowed_X_fft(windowed_X_fft==0) = 1e-5;
windowed_X_fft_dB_mag=20*log10(abs(windowed_X_fft));
% Rescale so that max. val. is 96 dB
windowed_X_fft_dB_mag = windowed_X_fft_dB_mag + (96-max(windowed_X_fft_dB_mag));

X=windowed_X_fft_dB_mag;

end