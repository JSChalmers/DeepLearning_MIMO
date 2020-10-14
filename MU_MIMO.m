clc
clear all
close all



num_subframe = 500000
num_sym = 100000

SNR = 0:2:27
for i = 1:length(SNR)
    if SNR(i)>19
        num_sym=10000;
    end
    snr_linear = 10^(SNR(i)/10);
    sigma = sqrt(0.5/snr_linear);
    
    error = 0;
    parfor frame = 1: num_subframe
    
        sym = randi([0,3], 2, num_sym);
        x_sym = qammod(sym, 4,'UnitAveragePower', true);
        
        H = (rand(2,2)+1j * randn(2,2))/sqrt(2);
       
        x_tr = inv(H)*x_sym;
        
        
        Es = mean(mean(abs(x_tr).^2));
        
        x_tr_norm = x_tr/sqrt(Es);
        
        
        noise = sigma*(randn(2,num_sym)+1j*(randn(2,num_sym)));
        x_rc = H * x_tr_norm + noise;
        
        
       
        
        to_decoder = x_rc * sqrt(Es);%% Normalize, such that can be decoded
        
        de_sym = qamdemod(to_decoder, 4, 'UnitAveragePower', true);
        
        correct = sum(sum((sym==de_sym),2));
        
        err = 2*num_sym - correct;
        error = error+err;
           
    end
    
    SER(i) = error/(2*num_sym*num_subframe)
    
end

figure()
semilogy(SNR, SER,'-*')
grid on
xlabel('SNR')
ylabel('SER')