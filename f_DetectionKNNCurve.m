
function s = f_DetectionKNNCurve(v_coeff)

%%%%% last update 31/07/2019

% vérifier les résultats avec normalisation

HOSth={

'   '  'B-PSK' 'Q-PSK' '8-PSK' '4-ASK' '8-ASK' '16-QAM'    ;
'M40'    1      -1       0       1.64    1.77    -0.67     ;
'M41'    1       0       0       1.64    1.77     0        ;
'M42'    1       1       1       1.64    1.77     1.32     ;
'M60'    1       0       0       2.92    3.62     0        ;
'M61'    1      -1       0       2.92    3.62    -1.32     ;
'M62'    1       0       0       2.92    3.62     0        ;
'M63'    1       1       1       2.92    3.62     1.96     ;
'M84'    1       1       1       5.25    7.92     3.12     ;
'C40'   -2      -1       0      -1.36   -1.42    -0.68     ;
'C41'   -2       0       0      -1.36   -1.42     0        ;
'C42'   -2      -1      -1      -1.36   -1.42    -0.68     ;
'C60'    16      0       0       8.32    7.19     0        ;
'C61'    16      4       0       8.32    7.19     2.08     ;
'C62'    16      0       0       8.32    7.19     0        ;
'C63'    16      4       4       8.32    7.19     2.08     ;

};

HOSth=cell2mat(HOSth(2:end,2:end));

lSNR=length(SNR);
lM=length(M);
lMMC=lM*lMC;

stat=zeros(4+Nt,lSNR*lMMC);

count=1;
for cptSNR=1:lSNR
    for cptM=1:lM
        for cptMC=1:lMC           
            stat(1,count)=SNR(cptSNR);
            stat(2,count)=M(cptM);
            count=count+1;      
        end
    end
end

% HOSth=HOSth(9:end,:);        %%% uniquement cum

for cptNt=1:Nt

% mat_HOS_cptNt=mat_HOS(15:end,:,cptNt); % HOC
% mat_HOS_cptNt=mat_HOS([1 3 5 6 8 10 12 14 15 16 17 18 19 20 21],:,cptNt); %
% HOS
mat_HOS_cptNt=mat_HOS([1 2 4 6 7 9 11 13 15 16 17 18 19 20 21],:,cptNt); %
% DEN

[~,resUz]=min(pdist2(HOSth',mat_HOS_cptNt','euclidean'));
stat(2+cptNt,:)=resUz;

end

if Nt==1
stat(4,:)=stat(3,:);
else
% stat(3+Nt,:)=mode(stat(3:2+Nt,:)); %nbr d'occurences
stat(3+Nt,:)=round(sum(v_coeff.*stat(3:2+Nt,:)));
end
stat(4+Nt,:)=double(stat(2,:)==stat(3+Nt,:));

courbe=zeros(lSNR,2);
courbe(:,1)=SNR;
end_=0;
for cptSNR=1:lSNR
    beg_=end_+1;
    end_=beg_+lMMC-1;
    courbe(cptSNR,2)=sum(stat(4+Nt,beg_:end_))/lMMC;
end

s = sum(courbe(:,2));

% figure
% plot(courbe(:,1),courbe(:,2))
