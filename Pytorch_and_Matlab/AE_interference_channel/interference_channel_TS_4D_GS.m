clc
close all
clear all
fid = fopen(['w4_256.txt']);
data_read_cell = textscan(fid , '', 'Delimiter', '\t', 'EmptyValue', NaN);
fclose(fid);
const = cell2mat(data_read_cell);
size(const)
const =sqrt(2)* const/sqrt(mean(sum(const.^2,2)));
mean(sum(const.^2,2))

M = 16^2;

SNR=0:1:14; %SNR in dB, average received power at one Rx element over the average noise power at that element
R = 4
for k=1:length(SNR)
    SNR(k)
    minErr=2e3;
    maxSym = minErr * 1e4;
    minSym = minErr * 1e2;
    totErr=0;
    totSym=0;
  
    while totErr<minErr & totSym<maxSym | totSym<minSym
        A=randi(M, 1,1); %transmitted alphabet
    
        st = const(A,:);
    
        snr_lin = 10^(SNR(k)/10); %dB to linear scale
        sigma = 0.5/(R * snr_lin) ; %the sigma square of the noise
        noise = sqrt(sigma)*randn(size(st)); %noise matrix 
        
        
        %Transceiver
        sr = st + noise;
          
        sr_rep = repmat(sr, M,1);
    
        

        
        [~, m_hat] = min(sum(abs(sr_rep-const).^2, 2));
  
        
        totErr = totErr + 1 - sum(eq(A,m_hat));
        totSym = totSym + 1;
        
    end
%     
    ser(:, k) =totErr/totSym;
%     
end



figure()
semilogy(SNR,ser, '-*', 'DisPlayName','Time-sharing, W4_256')
xlabel('SNR')
ylabel('BLER')
grid on
legend()
