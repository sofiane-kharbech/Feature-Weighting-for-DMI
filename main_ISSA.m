
clear
close all
clc

load base_adcom

%% Tested algorithms (Note : pos_xxx contains best values of x (i.e., C and sigma))

tic

% [mbest_SSA,stdbest_SSA,sem_SSA,mFEs_SSA,SR_SSA,pos_SSA,cg_curve_SSA] = f_SSA(40,500,30,M,mat_HOS,SNR,Nt,lMC); % FEs = npop*max_iter

[mbest_ISSA,stdbest_ISSA,sem_ISSA,mFEs_ISSA,SR_ISSA,pos_ISSA,cg_curve_ISSA] = f_ISSA(40,500,30,M,mat_HOS,SNR,Nt,lMC); % FEs = npop*max_iter

% [mbest_WSSA,stdbest_WSSA,sem_WSSA,mFEs_WSSA,SR_WSSA,pos_WSSA,cg_curve_WSSA] = f_WSSA(40,500,30,M,mat_HOS,SNR,Nt,lMC); % FEs = npop*max_iter

% [mbest_STSSA,stdbest_STSSA,sem_STSSA,mFEs_STSSA,SR_STSSA,pos_STSSA,cg_curve_STSSA] = f_STSSA(40,500,30,M,mat_HOS,SNR,Nt,lMC); % FEs =

toc