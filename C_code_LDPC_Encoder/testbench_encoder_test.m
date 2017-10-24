clear all;
mex LDPCEncoder_test.c;

B=3840; %input sequence length
c=int32(randi([0 1],1,B));     %info bits
BG=2;   %base graph
rate=1/5;

% Table of possible lifting sizes
Lift_size=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 20 22 24 26 28 30 32 36 ...
    40 44 48 52 56 60 64 72 80 88 96 104 112 120 128 144 160 176 192 ...
    208 224 240 256 288 320 352 384];

% Determine number of bits in codeword
if BG == 1
    Kb = 22;  % infomation bits
elseif BG == 2
    if B > 640
        Kb = 10;
    elseif B > 560
        Kb = 9;
    elseif B > 192
        Kb = 8;
    else
        Kb = 6;
    end
end

% find the minimum value of in all sets of lifting sizes in Table 
% denoted as Zc, such that Kb*Zc >= B,
Zc = min(Lift_size(Lift_size >= (B/Kb)));
if BG == 1
    padding_size = (22*Zc-B); % zero padding
    c_padding=[c zeros(1,padding_size)]; % bit sequence with zero padding
elseif BG==2
    padding_size = (Kb*Zc-B); % zero padding
    c_padding=[c zeros(1,padding_size)]; % bit sequence with zero padding
end

if BG==1
     if Zc==2||Zc==4||Zc==8||Zc==16||Zc==32||Zc==64||Zc==128||Zc==256
        load('no_shift_values_BG1_a=2.mat')
        load('Gen_shift_values_BG1_a=2.mat')
        load('pointer_shift_values_BG1_a=2.mat')
     elseif Zc==13||Zc==26||Zc==52||Zc==104||Zc==208
        load('no_shift_values_BG1_a=13.mat')
        load('Gen_shift_values_BG1_a=13.mat')
        load('pointer_shift_values_BG1_a=13.mat')
     end
     
elseif BG==2
    if Zc==2||Zc==4||Zc==8||Zc==16||Zc==32||Zc==64||Zc==128||Zc==256
        load('no_shift_values_BG2_a=2.mat')
        load('Gen_shift_values_BG2_a=2.mat')
        load('pointer_shift_values_BG2_a=2.mat')
    elseif Zc==3||Zc==6||Zc==12||Zc==24||Zc==48||Zc==96||Zc==192||Zc==384
        load('no_shift_values_BG2_a=3.mat')
        load('Gen_shift_values_BG2_a=3.mat')
        load('pointer_shift_values_BG2_a=3.mat')
     elseif Zc==7||Zc==14||Zc==28||Zc==56||Zc==112||Zc==224
        load('no_shift_values_BG2_a=7.mat')
        load('Gen_shift_values_BG2_a=7.mat')
        load('pointer_shift_values_BG2_a=7.mat')
    end
end

%tic
v=LDPCEncoder_test(c_padding, BG, Zc, No, R, p, Kb);
%toc
load('Gen_BG2_Z384_Kb10.mat')
v_t=mod(double(c_padding)*G,2);
isequal(v,v_t)
