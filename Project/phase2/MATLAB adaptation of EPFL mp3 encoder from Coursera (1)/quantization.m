function x = quantization(sample, sf, ba, QCa, QCb)
    % The main script will handle the determination of the required 
    % number of bits and other parameters for each subband.
    
    % This function will take five input arguments and it will return the quantized sample
    % uniformly quantized sample
    x=((QCa*sample/sf)+QCb)*2.^(ba-1);
end


