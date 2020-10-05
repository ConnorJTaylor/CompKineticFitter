function [TotalCombinationsFunctions,NumberofKValues] = GenerateKineticLaws_rewritten(KeepV,KeepVR,KeepVP,a,TotalCombinations,Combinations,Conc)
format short
Counter = 1;
Alphabet = ['A';'B';'C';'D';'E';'F';'G';'H'];
OrderedRateLaws = [];
OrderedRateLaw = []; %used as a temporary variable.
TotalCombinationsFunctions1 = {'1','2','3','4','5','6','7','8'};
TotalCombinationsFunctions = TotalCombinationsFunctions1(:,1:size(Conc,2));
Concentrations = [1,2,3,4,5,6,7,8];
Concentrations = Concentrations(1:size(Conc,2));
NumberofKValues = [];
%TotalFunctions = cell(TotalCombinations,size(Conc,2));
%changed from @(KV,a)0
startoffunction = '@(KV,a)0';
TotalFunctions = strings(TotalCombinations,size(Conc,2));
TotalFunctions(:,:) = startoffunction;

for i = 1:size(KeepV,2) %for each Combinations.CombX (X = 1,2,etc.).
    if Counter >= TotalCombinations
        break %stop looking for Combinations.CombX where X is greater than the KeepV size
    end 
    eval(strcat(sprintf('CombinationSize = size(Combinations.Comb%d,1);',i)))
    
    for L = 1:CombinationSize %for each combination of KeepV columns.
        clear KValues
        KValues = struct('Placeholder',1);
        eval(strcat(sprintf('CombColumns = Combinations.Comb%d(L,:);',i)))
        
        NewKeepV = [];
        NewKeepVR = [];
        NewKeepVP = [];
        for X = 1:size(CombColumns,2)
            NewKeepV(:,X) = KeepV(:,CombColumns(X)); %only the relevant reactions from KeepV.
            NewKeepVR(:,X) = KeepVR(:,CombColumns(X));
            NewKeepVP(:,X) = KeepVP(:,CombColumns(X));
        end
        
        
        
        for j = 1:size(NewKeepV,2) %for each of the columns, each has a specific k value (j)
            %work out the base reaction for the k value
            CurrentFunction = sprintf('KV(%d)',j);
            ReactantPositions = find(NewKeepVR(:,j)>0);
            ReactantPositions = ReactantPositions';
            
            for jj = 1:size(ReactantPositions,2)
                %this is the base reaction for this column and this k value
                %removed a.%s
            CurrentFunction = [CurrentFunction,sprintf('*a.%s^',Alphabet(ReactantPositions(jj))),num2str(NewKeepVR(ReactantPositions(jj),j))];
            end
            %if it's zero order, replace 9 with zero
            CurrentFunction = strrep(CurrentFunction,'9','0');
            
            for K = 1:size(Conc,2) %look at how each of the species are affected, species (K)
                %reactants
                if NewKeepVR(K,j) >0 && NewKeepVR(K,j) ~= 9 %if it reacts and is not 9
                    CurrentTotalFunction = char(TotalFunctions(Counter,K)); %must convert to character array first
                    CurrentTotalFunction = [CurrentTotalFunction,'-',num2str(NewKeepVR(K,j)),'*',CurrentFunction]; %then add characters
                    TotalFunctions(Counter,K) = CurrentTotalFunction; %then convert back to string
                elseif NewKeepVR(K,j) == 9 %if it is 9 (zero order)
                    CurrentTotalFunction = char(TotalFunctions(Counter,K)); %must convert to character array first
                    CurrentTotalFunction = [CurrentTotalFunction,'-',CurrentFunction]; %then add characters
                    TotalFunctions(Counter,K) = CurrentTotalFunction; %then convert back to string
                end
                
                %products
                if NewKeepVP(K,j) >0
                    CurrentTotalFunction = char(TotalFunctions(Counter,K)); %must convert to character array first
                    CurrentTotalFunction = [CurrentTotalFunction,'+',num2str(NewKeepVP(K,j)),'*',CurrentFunction]; %then add characters
                    TotalFunctions(Counter,K) = CurrentTotalFunction; %then convert back to string
                end
            end
        end
        NumberofKValues(Counter) = size(CombColumns,2);
        sprintf('Number of functions added: %d of %d',Counter,TotalCombinations)
        Counter = Counter+1;
        
    end
end


%turn to functions
TotalCombinationsFunctions = cellstr(TotalFunctions);

end