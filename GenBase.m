
%%%%% last update 07/11/2019

% clk=clock;
% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',str2num(strrep(num2str(clk(4:6)),' ',''))*1e3))

% ! sudo chmod -R 777 /home/sofiane/.matlab/local_cluster_jobs/R2017a
% parpool('local',3)

rng shuffle

clc
clear

tic

%%%%%%%%%%%%%%%%%%%%%%%

Nt = 2;
Nr = 6;
K = 2000;

M = [1 2 3 4 5 6];

% SNR = -10:5:20;
SNR = 2; %-10:1:15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lMC = 2000;

eestd = 0;
cfo = 0; % 1e-4;
phznoise = 0; % comm.PhaseNoise('Level',-3,'FrequencyOffset',2e-3,'SampleRate',1);

mat_HOS=[];

for cptSNR=1:length(SNR) %------------------------> Boucle SNR

    for cptM=1:length(M) %------------------------> Boucle sur la modulation
        
        Pb_ee = eestd*randn(1,lMC);
        parfor cptMC=1:lMC %------------------------> Boucle MC
                
            HOS=f_HOS_Extraction(M(cptM),SNR(cptSNR),K,Nt,Nr,Pb_ee(cptMC),cfo,phznoise);
            mat_HOS=[mat_HOS; permute(HOS,[3 2 1])];
                
        end
            
        % clc
        % pourc=(((cptSNR-1)*length(M)+cptM)/(length(M)*length(SNR)))*100;
        % disp(['Input Processing at ---> ' num2str(pourc) ' %'])
        disp(['M = ' num2str(M(cptM)) ' completed.'])
            
    end
end

mat_HOS=permute(mat_HOS,[2 1 3]);

clearvars -except M SNR Nt Nr lMC mat_HOS

toc