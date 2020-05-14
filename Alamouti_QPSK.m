clc
clear
close all
% N=500000 %total number of symbol pairs to be transmitted (should be at least 10 times more than expected 1/min(BER))
M=4; %PSK order (must be a power of 2): 2, 4, 8 etc'
Tx=2; %number of Tx elements, must be 2
Rx=1; %number of Rx elements
SNR=0:2:22; %SNR in dB, average received power at one Rx element over the average noise power at that element


for k=1:length(SNR)
    SNR(k)
    N=50000; %% Use more samples to get a smoother curve
    A=floor(M*rand(2,N)); %transmitted alphabet
    st = qammod(A,M, 'UnitAveragePower', true);
 
    snr=10^(SNR(k)/10); %dB to linear scale
    
    sig1=0.5/snr ; %the sigma square of the noise
    
    Ns=sqrt(sig1)*(randn(2*Rx,N)+j*randn(2*Rx,N)); %noise matrix 
    
    %Transceiver
    for n=1:N
        H=[]; %equivalent channel matrix initialization
        for r=1:Rx
            h=(randn(1,2)+j*randn(1,2))/sqrt(2); %Rayleigh channel
            H=[H; h(1) h(2); h(2)' -h(1)'];
        end %m
    
        sr(:,n)=H'*H*st(:,n)+H'*Ns(:,n); %received symbols
        temp = H'*H;
        sr(:,n) = sr(:,n)/temp(1,1);
    end %n
    
    decoded = qamdemod(sr, M,'UnitAveragePower',true );
    temp = xor(A-decoded, 0);
    temp = sum(temp);
    temp(temp>1)=1;
    ser(:, k) =sum(temp)/N;
        
end %k (SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(0:2:22,ser)
