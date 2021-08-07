
function [mbest,stdbest,sem,mFEs,SR,pos,cg_curve] = f_STSSA(N,Max_iter,RUNS,M,mat_HOS,SNR,Nt,lMC)

FEs = 0;
runs = 1;
lb = zeros(15,1);
ub = 1e3*ones(15,1);
dim = length(lb);
accep_level = .08;
CostFunction = @(x,M,mat_HOS,SNR,Nt,lMC) f_ObjFct(x,M,mat_HOS,SNR,Nt,lMC); % Ojbective function


while runs <= RUNS
%Initialize the positions of salps
SalpPositions=initialization(N,dim,ub,lb);


FoodPosition=zeros(1,dim);
FoodFitness=inf;


%calculate the fitness of initial salps

for i=1:size(SalpPositions,1)
    SalpFitness(1,i)=CostFunction(SalpPositions(i,:),M,mat_HOS,SNR,Nt,lMC);
    FEs = FEs + 1;
    all_cost(FEs) = SalpFitness(1,i);
    cg_curve(runs,FEs) = SalpFitness(1,i); 
end

[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);

for newindex=1:N
    Sorted_salps(newindex,:)=SalpPositions(sorted_indexes(newindex),:);
end

FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while l<Max_iter+1
    
    c1 = 2*exp(-(4*l/Max_iter)^2); % Eq. (3.2) in the paper
   
    for i=1:size(SalpPositions,1)
      
        SalpPositions= SalpPositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            
            SalpPositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        SalpPositions= SalpPositions';
    end
   
    for i=1:size(SalpPositions,1)
        
        Tp=SalpPositions(i,:)>ub';Tm=SalpPositions(i,:)<lb';SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
       
        r=rand();
        SalpPositions(i,:) = r*(ub'+lb')-SalpPositions(i,:); % % Eq. (3.4) in the paper 
        
        % Check if solutions go outside the search spaceand bring them back
        Tp=SalpPositions(i,:)>ub';Tm=SalpPositions(i,:)<lb';SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
          
        SalpFitness(1,i)=CostFunction(SalpPositions(i,:),M,mat_HOS,SNR,Nt,lMC);
        FEs = FEs + 1;
        all_cost(FEs) = SalpFitness(1,i);
     
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
        end
        
        cg_curve(runs,FEs) = FoodFitness; 
    end
     
   l = l + 1;
end

best(runs) = FoodFitness;
pos(runs,:) = FoodPosition;

if find(all_cost<=accep_level)
    fes(runs) = min(find(all_cost<=accep_level));
else
  fes(runs) = 0;
end

disp(['This is run number ' num2str(runs)]);

FEs = 0;
runs = runs+1;

end % end runs

% bbest = min(best);
mbest = mean(best);
% wbest = max(best);
stdbest = std(best);
sem = stdbest/sqrt(RUNS);

if fes==0
    mFEs = -1;
else
    indx = find(fes==0);
    fes(indx) = [];
    mFEs = mean(fes);
end

SR = length(find(fes~=0))/RUNS*100;
if RUNS >1
    cg_curve = mean(cg_curve);
end


end % end STSSA

%% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,1); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end
