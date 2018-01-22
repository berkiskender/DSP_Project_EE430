

close all
clc;

% Test Sample Length 
N=1024;

% Test Sample Vector 
xx=1:N;

% Signal Vector
sig=sin( ( pi/64)*xx  );
sig=sig';
sig_res=0*sig;


% MDCT window length
mdct_len=128;


% TOTAL Frame frame length
total_fr=2*(length(xx)/mdct_len )-1;

m_fr=zeros(mdct_len/2,total_fr);

mlh=mdct_len/2;
xf=1:mdct_len;
xr=1:mdct_len/2;


%%% Since Matlab starts index at 1, we have to use (xf-0.5) instead
%%% of (xf+0.5). Please refer to Wikipedia MDCT.  
wf=sin( (xf-0.5)*pi/ (mdct_len));



% MDCT process with 50% overlap
for r=1:total_fr
    
    m_fr(:,r)=mdct4(   wf(:).*sig( ((mdct_len/2)*(r-1)+1 ) : ((mdct_len/2)*(r-1)+mdct_len ) )    );
    
end

% IMDCT process with 50% overlap
for r=1:total_fr
    
   sig_res( ((mdct_len/2)*(r-1)+1 ) : ((mdct_len/2)*(r-1)+mdct_len ) )=   sig_res( ((mdct_len/2)*(r-1)+1 ) : ((mdct_len/2)*(r-1)+mdct_len ) )   +  wf(:).* imdct4( m_fr(:,r) );
   
end

% DIfference Calculation
diff = sig_res- sig;
d_len = length(diff);

figure
imagesc(m_fr)
title(' TF Domain of MDCT frames ') 

figure;
subplot(311)
plot(xx, sig);
title('input SIGNAL : SINE ')
axis ([ 1 N  -1 1])
subplot(312)
plot(xx, sig_res);
title('Restored SIGNAL ')
axis ([ 1 N  -1 1])
axis tight
subplot(313)
plot(1:d_len, diff);
axis ([ 1 N  -1 1])
title('Difference (Only have a difference in the non-overlapped perod) ')


figure;
plot(xf, wf);
title('Sine Window')




