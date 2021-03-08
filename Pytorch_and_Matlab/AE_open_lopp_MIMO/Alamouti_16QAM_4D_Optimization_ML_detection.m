clc 
close all
clear all

fid = fopen(['ae4_256.txt']);
data_read_cell = textscan(fid , '', 'Delimiter', '\t', 'EmptyValue', NaN);
fclose(fid);
const = cell2mat(data_read_cell); 

const_t1 = zeros(length(const), 2);
const_t1(:,1) = const(:,1) + 1j* const(:,3);  %%% create s1_t1
const_t1(:,2) = const(:,2) + 1j* const(:,4);  %%% create s2_t1

const_t2 = zeros(length(const), 2);
const_t2(:,1) = -const(:,2) + 1j* const(:,4); %%% create s1_t2
const_t2(:,2) = const(:,1) - 1j* const(:,3); %%% create s2_t2

const_t1 = const_t1.'; %%% Transpose to 2*256
const_t2 = const_t2.'; %%% Transpose to 2*256


M=16; %
Tx=2; %number of Tx elements, must be 2
Rx=1; %number of Rx elements
SNR=0:2:22; %SNR in dB, average received power at one Rx element over the average noise power at that element

for k=1:length(SNR)
    minErr=2e3;
    maxSym = minErr * 1e4;
    minSym = minErr * 1e2;
    totErr=0;
    totSym=0;
  
    while totErr<minErr & totSym<maxSym | totSym<minSym
        A=randi(M^2, 1,1); %transmitted alphabet
    
        st_t1 = const_t1(:,[A]);
        st_t2 = const_t2(:,[A]);
    
        snr=10^(SNR(k)/10); %dB to linear scale
        sigma=0.5/snr ; %the sigma square of the noise
        Ns=sqrt(sigma)*(randn(2*Rx,1)+j*randn(2*Rx,1)); %noise matrix 
        
        %Transceiver
        h=(randn(1,2)+j*randn(1,2))/sqrt(2); %Rayleigh channel
        sr_t1(:,1)=h*st_t1(:,1)+Ns(1,1); %received symbols
        sr_t2(:,1)=h*st_t2(:,1)+Ns(2,1); %received symbols
          
        sr_rep = repmat([sr_t1;sr_t2], 1, M^2); 
        
        ref_t1 = const_t1(:, [1:1:M^2]);
        ref_t2 = const_t2(:, [1:1:M^2]);
        
        ref = [h*ref_t1; h*ref_t2];
        
        [~, m_hat] = min(sum(abs(sr_rep-ref).^2, 1));
        
        totErr = totErr + 1 - sum(eq(A,m_hat));
        totSym = totSym + 1;
        
    end
%     
    ser(:, k) =totErr/totSym;
%     
end



figure()
semilogy(0:2:22,ser, '-*', 'DisPlayName','Alamouti w/ 4D GS ML detection')
xlabel('SNR')
ylabel('BLER')
grid on
legend()
