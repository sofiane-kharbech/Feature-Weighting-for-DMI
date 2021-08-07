
%%%%% last update 07/11/2019

function HOS=f_HOS_Extraction(M,SNR,K,Nt,Nr,pbee,cfo,phznoise)

% warning('off','comm:system:warnobsolete:obsoleteReplace'); % [~, id]= lastwarn
% warning('off','comm:obsolete:randint');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Emetteur

switch M
   case 1
      x = randi([0,1],1,Nt*K);
      x_mod = pskmod(x,2);
   case 2
      x = randi([0,3],1,Nt*K);
      x_mod = pskmod(x,4);
   case 3
      x = randi([0,7],1,Nt*K);
      x_mod = pskmod(x,8);
   case 4
      x = randi([0,3],1,Nt*K);
      x_mod = pammod(x,4);
   case 5
      x = randi([0,7],1,Nt*K);
      x_mod = pammod(x,8);
   case 6
      x = randi([0,15],1,Nt*K);
      x_mod = qammod(x,16);
   otherwise
      disp('Unknown modulation type');
end

scale = modnorm(x_mod,'avpow',1);   % Facteur de normalisation
x_mod = x_mod*scale;                % Normalisation de la puissance

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Canal

H=randn(Nr,Nt)+1i*randn(Nr,Nt);
y_mimo=H*reshape(x_mod,Nt,[]);

cfo_mat = repmat(exp(2i*pi*(1:1:K)*cfo),Nr,1);
y_mimo = y_mimo.*cfo_mat;
if phznoise ~=0
   for cptNr=1:Nr
       y_mimo(cptNr,:) = phznoise(y_mimo(cptNr,:).').';
   end
end

y_mimo_SNR = awgn(y_mimo,SNR,'measured');

Pb = mean(abs(reshape(y_mimo_SNR-y_mimo,1,[]).^2)) + pbee;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCMA

R = mean(real(x_mod).^4) / mean(real(x_mod).^2);
[W,yeg] = f_SCMA(Nt,Nr,y_mimo_SNR,R);
Pbf=Pb*W.'*conj(W);

% bf=W.'*(y_mimo_SNR-y_mimo)    % pour vérification
% mean(abs(bf(1,:).^2))         % pour vérification

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOS

for cpt=1:Nt
    HOS(cpt,:)=f_CalcHOS(yeg(cpt,:),Pbf(cpt,cpt));
end
