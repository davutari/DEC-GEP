%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

% CostFunction: Cost Function

% nVar: Number of Decision Variables



%VarMin: Lower Bound of Decision Variables
% VarMax: Upper Bound of Decision Variables

%% DE Parameters

% MaxIt: Maximum Number of Iterations

% nPop: Population Size

% beta_min: Lower Bound of Scaling Factor
% beta_max: Upper Bound of Scaling Factor

% pCR=0.2: Crossover Probability
% MesFlag: 1 for turning on the message

function [SolBest, SolCost, BestCost]=de(CostFunction,nVar,VarMin,VarMax,MaxIt,nPop,beta_min,beta_max,pCR,MesFlag)

VarSize=[1 nVar];   % Decision Variables Matrix Size
empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop    
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position,1,MaxIt);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
		y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position,it,MaxIt);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    if MesFlag==1
    disp(['Iteration: ' num2str(it) ': Best Cost: = ' num2str(BestCost(it))]);
    end    
end
SolBest=BestSol.Position;
SolCost=BestSol.Cost;