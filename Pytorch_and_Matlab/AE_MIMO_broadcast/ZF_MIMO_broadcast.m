clc
clear all
close all




M=4;
SNR = 1:4:40
SER = zeros(size(SNR));


for i = 1:length(SNR)
    if SNR(i)>19
        num_sym=10000;
    end
    
    snr = SNR(i);
    snr_linear = 10^(snr/10);
    sigma = sqrt(0.5/snr_linear);
    

    minErr=5e6;
    maxSym=minErr * 1e5;
   
    maxH=1e6;
    totErr=0;
    totSym=0;
    totH=0;
    num_sym=1e4;
    while totErr<=minErr & totSym<=maxSym & totH<=maxH
        H = (randn(2,2)+1j * randn(2,2))/sqrt(2);
        
        
        mess0 = [repmat([0,1,2,3],1,4); reshape(repmat([0,1,2,3],4,1), 1,[])];
        s0 = qammod(mess0, M,'UnitAveragePower', true);
        x_p = inv(H)*s0;
        Es = mean(mean(abs(x_p).^2));
        
        m = randi([0,3], 2, num_sym);
        s = qammod(m, 4,'UnitAveragePower', true);
        x_tilde = inv(H)*s;
        
        tx = x_tilde/sqrt(Es);
        Es4 = mean(mean(abs(tx).^2));
        
        noise = sigma*(randn(2,num_sym)+1j*(randn(2,num_sym)));
        x_rc = H * tx + noise;
        
     
        to_decoder = x_rc * sqrt(Es);%% Normalize, such that can be decoded
        m_hat = qamdemod(to_decoder, 4, 'UnitAveragePower', true);
        
        correct = sum(sum((m==m_hat),2));
        
        
        totH = totH+1;
        totErr =totErr +2*num_sym -correct;
        totSym = totSym + 2*num_sym;
        
    end
    
    SER(i) = totErr/(totSym)
    
end

%a=openfig('results_tmp.fig');
%hold on
figure()
semilogy(SNR, SER,'-*','DisplayName', 'ZF-2')

grid on
xlabel('SNR')
ylabel('SER')