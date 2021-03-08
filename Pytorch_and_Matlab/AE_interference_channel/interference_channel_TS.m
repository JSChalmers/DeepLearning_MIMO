clc
clear all
close all

SEP_QAM2_16 = @(EsNo) 3*qfunc(sqrt(2*EsNo/10)) - 9/4*qfunc(sqrt(2*EsNo/10)).^2;

EbNo=4:1:14;
EbNo_new = EbNo;

EsNo_linear = 4 * 10.^(EbNo_new/10);
ser11 = SEP_QAM2_16(EsNo_linear);
ser = 1- (1-ser11).*(1-ser11)

semilogy(EbNo, ser)
grid on