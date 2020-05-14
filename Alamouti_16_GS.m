clc
clear
close all
% N=500000 %total number of symbol pairs to be transmitted (should be at least 10 times more than expected 1/min(BER))
M=16; %PSK order (must be a power of 2): 2, 4, 8 etc'
Tx=2; %number of Tx elements, must be 2
Rx=1; %number of Rx elements
SNR=0:2:22; %SNR in dB, average received power at one Rx element over the average noise power at that element





a = [1.0364679 ,  1.2304633,  -0.07451034, -0.59903705 , 0.00228362, -0.824625,-0.5506731 , -0.38115218, -1.1879416 ,  0.07522251,  0.5699033,   0.6780808, 0.0400181 ,  0.3104105,   0.70689005 ,-1.0836288];
b =[0.40785137, -0.32265544,  1.109166,   -0.9286522,   0.00257764,  0.9765041,-0.19888385 , 0.46027398, -0.46247852, -1.2538868 , -0.13857605 , 1.0809493,-0.58626485 , 0.4999902,  -0.85292757 , 0.28704563];

const = a + 1j*b;  %%% 16QAM with geometric shaping



for k=1:length(SNR)
    SNR(k)
    N=50000;
    A=floor(M*rand(2,N)); %transmitted alphabet
    st = genqammod(A,const);
    
    snr=10^(SNR(k)/10);
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
    
    decoded = genqamdemod(sr, const);
    temp = xor(A-decoded, 0);
    temp = sum(temp);
    temp(temp>1)=1;
    ser(:, k) =sum(temp)/N;
    

     
    
end %k (SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(0:2:22,ser)

