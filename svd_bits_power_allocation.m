clc
clear all
close all


num_subframe = 10000
num_sym = 1000
SNR = 0:2:27;

for temp1 = 1: length(SNR)
    error = 0;
    error1=0; % 16QAM
    error2=0; % QPSK + QPSK
    error3=0; % 8-QAM + BPSK
    
    if SNR(temp1)<6
        num_sym = 50;
    elseif SNR(temp1)<11
        num_sym = 2000
    else
        num_sym = 10000
    end
   
    for j = 1:num_subframe
        snr_linear = 10^(SNR(temp1)/10);
        sigma = sqrt(0.5/snr_linear);
        noise = sigma*(randn(2,num_sym)+1j*(randn(2,num_sym)));
        
        H = (rand(2,2)+1j * randn(2,2))/sqrt(2);
        [U,S,V] = svd(H);
        precoder=V;
        shaping=U';
       
   
        
        %%% Case1 : 16-QAM on one stream:
        Es=2;
        err_case1=0;
        sym = randi([0,15], 1, num_sym);
        tx_sym = sqrt(Es) * qammod(sym, 16, 'UnitAveragePower', true);
        rx_sym = shaping(1,:) * (H * precoder(:,1)*tx_sym + noise(1,:));

        rx_norm = rx_sym/sqrt(Es)/S(1,1);
    
        decoded = qamdemod(rx_norm,16,'UnitAveragePower', true);
        correct =sum((sym==decoded));
        err_case1 = num_sym - correct;
        error1 = error1 + err_case1;
       
        %%%% Case2: QPSK on both streams:
        err_case2=0;
        temp = 0.1:0.02:1;
        for index = 1:length(temp)
            es1 = temp(index);
            es2 = 2-es1;
           
            sym1 = randi([0,3], 1, num_sym);
            sym2 = randi([0,3], 1, num_sym);
            tx_sym1 = sqrt(es1) * qammod(sym1, 4, 'UnitAveragePower', true);
            tx_sym2 = sqrt(es2) * qammod(sym2, 4, 'UnitAveragePower', true);
        
            tx_sym = [tx_sym1;tx_sym2];
            rx_sym = shaping*(H*precoder* tx_sym + noise); 
            
            rx_sym1 = rx_sym(1,:);
            rx_norm1 = rx_sym1/sqrt(es1)/S(1,1);
            rx_sym2 = rx_sym(2,:);
            rx_norm2 = rx_sym2/sqrt(es2)/S(2,2);
           
            decoded1 = qamdemod(rx_norm1,4,'UnitAveragePower', true);
            decoded2 = qamdemod(rx_norm2,4,'UnitAveragePower', true);
            
            temp21=(sym1==decoded1);
            temp22=(sym2==decoded2);

            temp2 = temp21+temp22;
            corr2 = numel(find(temp2==2));
            err2 = num_sym-corr2;
            
            err_qpsk(index) = err2;
        end 
        err_case2 = min(err_qpsk);
        error2 = error2 + err_case2;
       
        
        
        %%%%%% Case3: 8QAM one stream, BPSK another stream %%%%
        sym1 = randi([0,7], 1, num_sym);
        sym2 = randi([0,1], 1, num_sym);

        const_8qam = [1+1j, -1+1j, -1-1j,1-1j, (sqrt(3)+1), (sqrt(3)+1)*1j, -(sqrt(3)+1), -(sqrt(3)+1)*1j];
        Es =  mean(abs(const_8qam).^2);
        const= const_8qam/sqrt(Es);
        
        err_case3=0;
        temp = 0.1:0.02:1;
        for index = 1:length(temp)
            es1 = temp(index);
            es2 = 2-es1;
            x_sym1 = sqrt(es1) * genqammod(sym1, const); 
            x_sym2 = sqrt(es2) * qammod(sym2, 2,'UnitAveragePower', true);
            
            x_tilde2 = [x_sym1;x_sym2];
            x_tr2 = precoder*x_tilde2;
            
            y2=H*x_tr2 + noise;
            y_tilde2=shaping*y2;

            decoded2 = inv(S)*y_tilde2;
            stream1 = decoded2(1,:)/sqrt(es1);
            stream2 = decoded2(2,:)/sqrt(es2);
            
            dec_sym1 = genqamdemod(stream1,const);
            dec_sym2 = qamdemod(stream2,2,'UnitAveragePower', true);

            temp21=(sym1==dec_sym1);
            temp22=(sym2==dec_sym2);

            temp2 = temp21+temp22;
            corr3 = numel(find(temp2==2));
            err3 = num_sym-corr3;
            err_8qam(index) = err3;  
        end
        err_case3 = min(err_8qam);
        error3 = error3+err_case3;
        
        err = min([err_case1, err_case2, err_case3]);
        error = error + err; 
    end
    
    SER(temp1) = error/(num_sym*num_subframe);
    SER_16QAM(temp1) = error1/(num_sym*num_subframe);
    SER_QPSK(temp1)= error2/(num_sym*num_subframe);
    SER3(temp1)= error3/(num_sym*num_subframe);
   
      
end
%%

figure()
semilogy(SNR, SER,'-o')

legend('bits-power-allocation')



