% Beijing Jiaotong University
% Author: Bai
% Contention Based Massive Access Protocol for B5G£º A Compressive Sensing Method
% Frequency-Selective Channel
% Channel bandwidth: 1.4MHz
% Subcarrier bandwidth: 15kHz
% Number of subcarriers: 6*12=72
% 128 FFT
% Protection of Bandwidth: 128-72=56

close all
clear all
clc
warning('off')

%====================System Parameter=========================
% No. of slot
NUM_slot = 3;
% No. of resource elements used for pilot and data transmission
m = 84*NUM_slot;
% Modulation:1-QPSK,2-16QAM,3-64QAM
Mod_Type = 3; 
%======================Drawing Figure=========================
% Matrix storing point for drawing figure
NUM_ROW  = 7; % the first row for num of act users
drawPointMtx1 = zeros(NUM_ROW,1);
%====================User Parameter====================
% No.of all users
NUM_ALL_USERS = 1000000;
% Length of Pilot in each RB*NUM_slot
mp = 72;
% No. of data symbol, each user in in each RB*NUM_slot
nd = 2; 
% Data code-word length
md = floor((m-mp)/nd); 

if Mod_Type==1 % QPSK symbols set
    A=[1+1j,1-1j,-1+1j,-1-1j];
    BIT_PER_SBL = 2;
elseif Mod_Type==2 % 16QAM symbols set
    A=[-3+3j,-1+3j,1+3j,3+3j,-3+1j,-1+1j,1+1j,3+1j,-3-1j,-1-1j,1-1j,3-1j,-3-3j,-1-3j,1-3j,3-3j];
    BIT_PER_SBL = 4;
elseif Mod_Type==3 % 64QAM symbols set
    A = [1+1j,1+3j,1+5j,1+7j,3+1j,3+3j,3+5j,3+7j,5+1j,5+3j,5+5j,5+7j,7+1j,7+3j,7+5j,7+7j,1-1j,1-3j,1-5j,1-7j,3-1j,3-3j,3-5j,3-7j,5-1j,5-3j,5-5j,5-7j,7-1j,7-3j,7-5j,7-7j,-1+1j,-1+3j,-1+5j,-1+7j,-3+1j,-3+3j,-3+5j,-3+7j,-5+1j,-5+3j,-5+5j,-5+7j,-7+1j,-7+3j,-7+5j,-7+7j,-1-1j,-1-3j,-1-5j,-1-7j,-3-1j,-3-3j,-3-5j,-3-7j,-5-1j,-5-3j,-5-5j,-5-7j,-7-1j,-7-3j,-7-5j,-7-7j];
    BIT_PER_SBL = 6;
end


drawPointMtxROW = 1; 
for NUM_PILOT = [64 200 500 1000 5000 10000 50000 ]    % No. of pilot [64 200 300 400 500] 
    drawPointMtxColumn = 2;
    drawPointMtxROW = drawPointMtxROW + 1;
    drawPointMtx1(drawPointMtxROW,1) = NUM_PILOT;
    for NUM_ACT_USERS = 5:5:100
        drawPointMtx1(1,drawPointMtxColumn) = NUM_ACT_USERS;   
        NUM_ITR = 10000;
        for itr = 1:NUM_ITR
            % Pilot code-word pool (non-orthometric)
            cwplt = A(randi(length(A),mp,NUM_PILOT));
            % Symbol code-word matrix (Corresponds to pilot)
            cwsbl = 2*randi([0,1],md,NUM_PILOT)-1;  
            %=================Step 1: random send pilot & ID &data  =================
            % Indices of selected pilot (may be repeated)
            pilotIndx = unidrnd(NUM_PILOT,[1,NUM_ACT_USERS]);
            % Indices of selected pilot (not repeated)        
            [pilot_choose, ia, ic] = unique(pilotIndx);
            %NO of collision pilot 
            pilot_collision = pilotIndx;
            pilot_collision(ia)=[];
            num_pilot_collision = length(unique(pilot_collision));            
            drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) + num_pilot_collision/NUM_ACT_USERS;              
        end
        drawPointMtxColumn = drawPointMtxColumn + 1;
    end
end

drawPointMtx1(2:NUM_ROW,2:21)  = drawPointMtx1(2:NUM_ROW,2:21) ./ NUM_ITR;
%save('figureprob','drawPointMtx1')
