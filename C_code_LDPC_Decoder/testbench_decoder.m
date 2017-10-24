clear all;
mex LDPCDecoder.c;
mex ML.c;
%mex test_generate_H_2.c;
no_iteration=int16(25);
%SNR_dB=-2.6:0.2:-2;
SNR_dB=30;
SNR=10.^(SNR_dB./10);
E=1;    %average energy of transmitted signal
muy=0;  %mean of noise
error=zeros(1,size(SNR_dB,2));
P=zeros(1,size(SNR_dB,2));

load('Shift_value_BG2_set1.mat')
load('Column_position_BG2.mat')
load('no_one_element_BG2.mat')

load('Gen_BG2_Z128_Kb10.mat')

Z=int16(128);
BG=int16(2);
rate=0.2;
Kb=int32(10);
no_punctured_columns=0;
rows=size(no_one_element,2)*Z;
cols=52*Z;

p_channel=zeros(1,cols);

for i2=1:size(SNR_dB,2)
    
std=sqrt(E/SNR(i2));    %std of noise

for m=1:10000
    u=randi([0 1],1,cols-rows);     %info bits
    v=mod(u*G,2);                  %codeword
    v_transmitted=zeros(1,size(v,2));
    v_transmitted(v==0)=sqrt(E/2);    %QPSK
    v_transmitted(v==1)=-sqrt(E/2);   
    w1=sqrt(E/(2*SNR(i2)))*randn(1,size((2*Z+1):(size(v,2)-no_punctured_columns*Z),2)); %noise in AWGN channel
    r=v_transmitted((2*Z+1):(size(v,2)-no_punctured_columns*Z))+w1;                  %received codeword 
    
%LLR of the channel
%        p_channel((2*Z+1):(size(v,2)-no_punctured_columns*Z))=(sqrt(2)/std^2).*r;             
%        p_channel_fixed=int16(round(p_channel*2^7));     %fixed point 9_7
    
   [p_channel]=ML(r,std,size(v,2),no_punctured_columns,Z);

   [v_estimate]=LDPCDecoder(no_iteration,shift_value,Col_position,no_one_element,Z,BG,Kb,rate,p_channel);

    if isequal(v,v_estimate)==0
        error(i2)=error(i2)+1;
    end
end

P(i2)=error(i2)/10000;

end
figure(1)
semilogy(SNR_dB,P);
xlabel('SNR_dB')
ylabel('BLER_dB');
