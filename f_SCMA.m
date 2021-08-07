
%%%%% last update 10/09/2017

function [G,yeg]=f_SCMA(Nt,Nr,y_mimo_SNR,R)

%--------------> Pr�-traitement de blanchiement (Batch) du signal re�u <--------------%

Ry = (y_mimo_SNR*y_mimo_SNR')/size(y_mimo_SNR,2);

[U,SIG] = eig(Ry);                          % D�composition en valeurs propres
SIG = SIG(Nr:-1:(Nr-Nt+1),Nr:-1:(Nr-Nt+1)); % Nr>=Nt
U = U(:,Nr:-1:(Nr-Nt+1));

B = (SIG^-0.5)*U'; % Matrice de blanchiement
YW = B*y_mimo_SNR; % Signal re�u blanchit

%--------------> Estimation du s�parateur W <--------------%        

W = eye(Nt);                       % Initialisation du s�parateur
%W = toeplitz(1:Nt);
%W = toeplitz(1:Nr,1:Nt);

mu=1e-2;                           % Pas de gradient stochastique
%mu=5e-3;

for ns=1:size(YW,2)                % Alg. de gradient stochastique

    YW_ns = YW(:,ns);
    z_ns = W.'*YW_ns;
    
    for n=1:Nt
 
        e = (real(z_ns(n)).^2-R).*real(z_ns(n));
        W(:,n) = W(:,n)-mu.*e.*conj(YW_ns);
        
    end
    
    %--------------> Post-traitement d'orthonormalisation <--------------%
    
    W=f_GramSchmidt(W); % Orthonormalisation � chaque it�ration ns
    
end


%% %%%%%%%%%% S�paration de source : Application du s�parateur W

yeg = W.'*YW;
G = B.'*W;