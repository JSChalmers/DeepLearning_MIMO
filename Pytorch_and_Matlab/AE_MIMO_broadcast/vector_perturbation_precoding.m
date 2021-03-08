clc
clear all
close all


M=4;
Nt=2;
Nr=2;

A=2;
modA  = @(x) x - A .* floor((x+A/2)/A);

SNR_dB = [1,5,9,13,17,21,25];%

SER = []
tic
for II = SNR_dB
    II
    SNR = 10^(II/10);
    sigma = sqrt(0.5/SNR);
    
    num_sym = 1e4
    Xi = 0.6;
    beta = Xi / SNR
    
    minErr=1e5;
    maxSym=minErr/1e-5;
    totErr=0;
    totSym=0;
    num_H=0;
    max_H = 1e4;
    
    while totErr<=minErr & totSym<=maxSym | num_ch < max_H
        num_H = num_H+1;
        H =  (randn(Nr, Nt) + 1j*randn(Nr, Nt)) / sqrt(2);
    
        x_tilde = H2symbol(H, beta, A);
        
        alpha = 1/ sqrt(mean(mean(abs(x_tilde).^2)));
        x_tilde = x_tilde *alpha;


        m = randi(2^M,1,num_sym);
        m1 = mod(m-1, M);
        m2 = floor((m-1)/M);
        Tx = x_tilde(:, m);

        Rx = H*Tx ;
        Rx = Rx + (randn(size(Rx)) + 1j*randn(size(Rx))) *sigma;

        rx1 = Rx(1,:)./alpha;
        rx2 = Rx(2,:)./alpha;



         %%%%% Receiver Operations
        mod_rx1 = modA(real(rx1)) + 1j*modA(imag(rx1));
        mod_rx2 = modA(real(rx2)) + 1j*modA(imag(rx2));

        ref = qammod([0,1,2,3], M)./2;

    %     figure()
    %     scatter(real(mod_rx2), imag(mod_rx2),20,'filled')
    %     hold on
    %     legend('mod-Rx')
    %     scatter(real(ref), imag(ref), 60,'filled')
    %     grid on


        rx1_rep = repmat(mod_rx1, [M,1]);
        rx2_rep = repmat(mod_rx2, [M,1]);

        ref = qammod([0,1,2,3], M)./2;
        ref_sym = repmat(ref.', [1,length(rx1_rep)]);

    %     figure()
    %     set(gcf,'position',[10,10,600,300])
    %     subplot(1,2,1)
    %     scatter(real(rx1), imag(mod_rx1))
    %     hold on
    %     scatter(real(ref), imag(ref), 40)
    %     subplot(1,2,2)
    %     scatter(real(rx2), imag(mod_rx2))
    %     hold on
    %     scatter(real(ref), imag(ref), 40)

        [~, dec1] = min(abs(rx1_rep - ref_sym).^2);
        ser1 = 1 - mean(eq(dec1-1, m1));

        [~, dec2] = min(abs(rx2_rep - ref_sym).^2);
         ser2 = 1 - mean(eq(dec2-1, m2));


        totErr =totErr +2*num_sym -sum(eq(dec1-1, m1)) -sum(eq(dec2-1, m2));
        totSym = totSym + 2*num_sym;
    %     ser = [ser_test1, ser_test2;ser1,ser2]
    end

   

    ser = totErr/totSym;
    SER = [SER, ser]
end
toc
 


openfig('results_tmp.fig')
hold on
semilogy(SNR_dB, SER,'-->')





function [out1, out2] = H2symbol(H,beta, A)
    Nt=2;
    Nr=2;
    M=4;
    W = H' * inv(H*H' + beta * eye(Nr));
    
    mess = [repmat([0,1,2,3],1,4); reshape(repmat([0,1,2,3],4,1), 1,[])];
    s = qammod(mess, M)./2;
    x_tilde = s;

    for I = 1:length(x_tilde)  
        
     %%%%%%    Transmitter
       
        
        Zr = -5:5;
        Zi = -5:5;
        [XX,YY] = ndgrid(Zr,Zi);
        Z = [XX(:),YY(:)];
        CZ = Z(:,1) + 1j*Z(:,2);

        p_star = A*[CZ(1);CZ(1)];
        Norm0 = sum(abs(W *(s(:,I) + p_star)).^2);
        for i=1:length(CZ)
            for j = 1:length(CZ)
                p_prime = A*[CZ(i); CZ(j)];
                coded = W *(s(:,I) + p_prime);
                Norm = sum(abs(coded).^2);
                if Norm <= Norm0
                    Norm0 = Norm;
                    p_star =p_prime;
                end
            end
        end
        
        x_tilde(:,I) = W *(s(:,I) + p_star);
    end
    
    alpha = sqrt(mean(mean((abs(x_tilde).^2))));
    out1 = x_tilde;
    out2 = alpha;

end


