# IWCMC_contention-based-CS-access-scheme
All matlab code of paper: Contention Based Massive Access Scheme for B5G: A Compressive Sensing Method


# Folder structure
figure_diff_acc: main code (different number of active users) that gennerate simulation results of Fig.3 and Fig.4 in paper
figure_diff_assump.m: main code (different number of active users) that gennerate simulation results of Fig.6 in paper
figure_diff_l.m: main code (different length of pilot) that gennerate simulation results of Fig.7 and Fig.8 in paper
figure_diff_snr.m: main code (different snr) that gennerate simulation results of Fig.9 in paper

LTE_FOUR_STEP: Function that performs pilot detection and channel estimation and data recovery in LTE RACH.
OMP_blind.m: Function that performs active user detection/pilot detection and channel estimation using OMP algorithm.
Data_LS.m: Function that performs data recovery including Demodulation、Deinterleaving、Decoding et.al.

figure_diff_acc_1000000user_l72.mat: Simulation results of Fig.3 and Fig.4 in paper
drawpicture3.m: Plot Fig.3 and Fig.4 in paper

prob.m: Simulations of Fig.5, generate the probability of collision under different number of active users and pilots
figureprob.mat: Simulation results of Fig.5 in paper

figure_difass.mat: Simulation results of Fig.6 in paper 
drawpicture6.m: Plot Fig.6 in paper 

figure_diff_snr_act35_l72_ue1000000.mat: Simulation results of Fig.7 and Fig.8 in paper  
drawpicture8.m: Plot Fig.7 and Fig.8 in paper 

figure_diff_l_1000000user_ac40.mat: Simulation results of Fig.9 in paper 
drawpicture9.m: Plot Fig.9 in paper 

# Citation
@INPROCEEDINGS{9148341,
  author={Bai, Yanna and Chen, Wei and Ai, Bo and Zhong, Zhangdui},
  booktitle={2020 International Wireless Communications and Mobile Computing (IWCMC)}, 
  title={Contention Based Massive Access Scheme for B5G: A Compressive Sensing Method}, 
  year={2020},
  volume={},
  number={},
  pages={1854-1859},
  doi={10.1109/IWCMC48107.2020.9148341}}
