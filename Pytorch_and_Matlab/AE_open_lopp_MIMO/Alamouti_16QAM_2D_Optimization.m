clc
clear
close all
% N=500000 %total number of symbol pairs to be transmitted (should be at least 10 times more than expected 1/min(BER))
M=16; %PSK order (must be a power of 2): 2, 4, 8 etc'
Tx=2; %number of Tx elements, must be 2
Rx=1; %number of Rx elements
SNR=0:2:22; %SNR in dB, average received power at one Rx element over the average noise power at that element


gs16 = [-0.7806694  -1.2768254   0.8599484  -0.08432804 -0.61991584  0.50889343 ...
   0.5126908   1.103399   -0.9956602  -0.03835376  0.11500008 -0.677602...
  -0.30583203  0.39765245  1.0609547  -0.16846679]...
+ 1j*[ 0.80480534  0.24740237 -0.9286406   0.5790462   0.12658182  0.95518804...
   0.23866254  0.5298605  -0.43862224 -0.04386776 -1.1040015  -1.1306624...
  -0.5748732  -0.45713484 -0.20805919  1.2538805 ];

gs16 = [6.7413801e-01 -3.7566015e-01  1.2275851e+00  1.0295780e+00...
   1.7410155e-02 -1.2162055e+00 -1.0585674e+00  5.9810090e-01...
   1.7279519e-02  3.5183391e-01 -5.6507850e-01 -7.7046734e-01...
  -6.2502789e-01  1.3913669e-02 -3.5433076e-02  7.1659893e-01]...
        +1j * [-8.7854159e-01  4.6747768e-01 -3.6506271e-01  3.6437711e-01...
  -1.2685099e+00 -4.1955590e-01  3.0192161e-01 -1.6512145e-01...
  -6.0401165e-01  5.1224715e-01 -2.0169437e-01  1.0017546e+00...
  -8.9360678e-01  7.8599673e-04  1.0975660e+00  1.0499746e+00];




for k=1:length(SNR)
    SNR(k)
    minErr=1e4;
    maxSym = minErr * 1e4;
    totErr=0;
    totSym=0;
   
    while totErr<minErr & totSym<maxSym | totSym<1e5
        N=1000; %% Use more samples to get a smoother curve
        A=floor(M*rand(2,N)); %transmitted alphabet
        st = genqammod(A,gs16);
        snr=10^(SNR(k)/10); %dB to linear scale
        sigma=0.5/snr ; %the sigma square of the noise
        Ns=sqrt(sigma)*(randn(2*Rx,N)+j*randn(2*Rx,N)); %noise matrix 
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

        decoded = genqamdemod(sr, gs16 );
        temp = xor(A-decoded, 0);
        temp = sum(temp);
        temp(temp>1)=1; 
        error =  sum(temp);
        totErr = totErr+error;
        totSym = totSym + N;
    end
    ser(:, k) =totErr/totSym;
        
end %k (SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(0:2:22,ser, '-*','DisplayName', 'Alamouti w/ 2D-GS-16QAM')
legend()
xlabel('SNR')
ylabel('BLER')
grid on