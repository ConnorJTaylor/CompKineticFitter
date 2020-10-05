function [ErrorTable,KTable,NoOfKTable] = OptimiseK_2(KeepV,KeepVR,KeepVP,Results2,Conc,Time,QuantifyWhich,catalytic)
% Main optimisation code.
% Version 2.3.5
% Last updated: 4th October 2020
% Connor Taylor
% University of Leeds

%%
%Initialisation
a = struct('A',0,'B',0,'C',0,'D',0,'E',0,'F',0,'G',0,'H',0); %alphabet struct.
clear Combinations %clear Combinations struct from any previous runs.
KV = []; %K Values definition
Combinations = struct('Placeholder',1); %a placeholder struct for all combinations of model terms.
TotalCombinations = 0; %early definition as a placeholder.
Counter = 0; %early definition as a placeholder.
delete datafile.mat %clear file from previous runs.
catalytic; %is the reaction a catalytic reaction

StartingMaterials = find(Conc>0); %used for rule 2, see README file for more information
%%
%find all rate laws
[AdjustedKeepVP,AdjustedKeepVR,NumberOfProdReactions] = AdjustKeepVforPowersNEW(KeepVP,KeepVR,catalytic);
AdjustedKeepV = AdjustedKeepVP-AdjustedKeepVR;


%%
%Rule 2 implementation
ContainsAStartingMaterial = [];
for i = 1:size(StartingMaterials,2)
    ContainsAStartingMaterial(i,:) = AdjustedKeepVR(StartingMaterials(i),:);
end

TempNecessaryReactionLogical = ones(1,size(ContainsAStartingMaterial,2));
for i = 1:size(ContainsAStartingMaterial,2)
    Currentterm1 = sum(ContainsAStartingMaterial(:,i));
    Currentterm2 = sum(AdjustedKeepVR(:,i))-Currentterm1;
    
    %if there are no starting materials in the model term, this reaction
    %cannot exist by itself. Therefore it's not a necessary reaction.
    if Currentterm1 < 1
        TempNecessaryReactionLogical(i) = 0;
    else
        %if the starting material is in the model term, but requires the
        %presence of another material that has not been generated yet, then
        %it's not able to exist by itself. Therefore it's not a necessary
        %reaction.
        if Currentterm2 >0
            TempNecessaryReactionLogical(i) = 0;
        end
    end
end


%%
%Compile reaction models & implement rules
for i = 1:4
    i
    %Every single possible combination of the columns are inputted into the
    %'Combinations' struct.
    CurrentCombinations = nchoosek(1:size(AdjustedKeepVR,2),i);
    KeepCombinationMatrix = zeros(size(CurrentCombinations,1),size(NumberOfProdReactions,1));
    KeepCombinationsLogical = ones(size(CurrentCombinations,1),1);
    
    for ii = 1:size(CurrentCombinations,1) %look through every combination
        TempComb = CurrentCombinations(ii,:); %specify the current combination of interest
        
        for iii = 1:size(KeepCombinationMatrix,2) %for each line of column combinations
            %find where the numbers are in the array that are between each
            %of the bounds
            NoOfOccurs = find(TempComb >= NumberOfProdReactions(iii,1) & TempComb <=NumberOfProdReactions(iii,2));
            %count the number of times there are numbers in each of the
            %ranges
            Number = size(NoOfOccurs,2);
            %assign the space in the matrix
            KeepCombinationMatrix(ii,iii) = Number;
            
            if Number > 1
                
                KeepCombinationsLogical(ii) = 0;
                break
            end
        end
        
        CurrentTotalNecessaryReactions = 0;
        for iii = 1:size(TempComb,2)
            CurrentTotalNecessaryReactions = CurrentTotalNecessaryReactions + TempNecessaryReactionLogical(TempComb(iii));
        end
        if CurrentTotalNecessaryReactions == 0
            KeepCombinationsLogical(ii) = 0;
        end
        
        if KeepCombinationsLogical(ii) == 1
            AvailableSpecies = find(Conc>0);
            UsedSpecies = [];
            %find products
            for iii = 1:size(TempComb,2)
                Tempavailable = find(AdjustedKeepVP(:,TempComb(iii))>0);
                Tempavailable = Tempavailable';
                AvailableSpecies = [AvailableSpecies Tempavailable];
                AvailableSpecies = unique(AvailableSpecies','rows');
                AvailableSpecies = AvailableSpecies';
            end
            
            AvailableSpeciesLogical = zeros(1,size(Conc,2));
            UsedSpeciesLogical = zeros(1,size(Conc,2));
            for iii = 1:size(AvailableSpecies,2)
                AvailableSpeciesLogical(AvailableSpecies(iii)) = 1;
            end
            
            %find reactants
            for iii = 1:size(TempComb,2)
                Tempused = find(AdjustedKeepVR(:,TempComb(iii))>0);
                Tempused = Tempused';
                UsedSpecies = [UsedSpecies Tempused];
                UsedSpecies = unique(UsedSpecies','rows');
                UsedSpecies = UsedSpecies';
            end
            
            for iii = 1:size(UsedSpecies,2)
                UsedSpeciesLogical(UsedSpecies(iii)) = 1;
            end
            
            %see if more species are used than available
            DifferenceLogical = AvailableSpeciesLogical-UsedSpeciesLogical;
            DifferenceLogicalMin = min(DifferenceLogical);
            
            if DifferenceLogicalMin < 0
                KeepCombinationsLogical(ii) = 0;
            end
        end
    end
    
    KeepCombinationsLogical = logical(KeepCombinationsLogical);
    CurrentCombinations = CurrentCombinations(KeepCombinationsLogical,:);
    
    eval(strcat(sprintf('Combinations.Comb%d',i),sprintf(' = CurrentCombinations;')))
    TotalCombinations = eval(sprintf('TotalCombinations + size(Combinations.Comb%d,1);',i));
end

%Function that compiles all possible rate laws from reaction models defined.
%Tip: Comment out the function and input 'TotalCombinationsFunctions' and
%'NumberofKValues' into 'OptimiseK_2' to avoid using this function.
[TotalCombinationsFunctions,NumberofKValues] = GenerateKineticLaws_rewritten(AdjustedKeepV,AdjustedKeepVR,AdjustedKeepVP,a,TotalCombinations,Combinations,Conc);

% try %try operator necessary in case parpool is already open.
%     core_info = eval('feature(''numcores'')');
%     parpool(core_info); %prepare for parfor loop. Replace with how many cores your computer has.
% catch
% end

%%
%Kinetic fitting stage

%Number of 'blocks' of 100 to be optimised.
NumberofFunctionArrays = ceil(TotalCombinations/100);

%Pre-allocate output tables for use in parfor.
%These can be defined by the use of 'NumberofFunctionArrays', such as:
%'ErrorTable = zeros(NumberofFunctionArrays,100)'
%or can be defined by however many optimisations you want to run, such as:
%'ErrorTable = zeros(150,100)'.
ErrorTable = zeros(NumberofFunctionArrays,100);
KTable = cell(NumberofFunctionArrays,100);
NoOfKTable = zeros(NumberofFunctionArrays,100);

interupt_time_fmincon = 200; %define maximum optimisation time.

%Define output function and fmincon options.
outputFunction = @(Time,Conc,state)interuptFunfmincon(Time,Conc,state,interupt_time_fmincon);
options = optimoptions('fmincon','Display','off','OutputFcn',outputFunction);


%Select the CombNo values you want to evaluate. For all evaluations use
%'CombNo = 1:NumberofFunctionArrays', or for example use 'CombNo = 23:24'
%to evaluate 2301-2500 within TotalCombinationsFunctions.
for CombNo = 1:NumberofFunctionArrays
    try
        %If there are less than 100 terms in the last block this will
        %error, which is why this is in a try operator.
        %CurrentCombinationsFunctions are the current functions within
        %the current block to be optimised.
        CurrentCombinationsFunctions = TotalCombinationsFunctions(((CombNo-1)*100)+1:((CombNo-1)*100)+100,:);
    catch
        %This code creates a 100 model block even if there aren't 100
        %models left.
        CurrentCombinationsFunctions = cell(100,size(TotalCombinationsFunctions,2));
        CurrentCombinationsFunctions(:,:) =({'@(KV,a)0'});
        Remaining = 100-((NumberofFunctionArrays*100)-size(TotalCombinationsFunctions,1));
        Bulk = size(TotalCombinationsFunctions,1)-Remaining;
        CurrentCombinationsFunctions(1:Remaining,:) = TotalCombinationsFunctions(Bulk+1:end,:);
    end
    
    
    for I = 1:100 %for every model within the block. Can be replaced with 'parfor'
        %The number of K values to be optimised is first established, from there
        %the bounds can be established as well as the starting point.
        try
            %This 'try' is necessary for when blocks have less than 100 terms in
            %them.
            CurrentNoofKValues = NumberofKValues((CombNo-1)*100+I);
        catch
            CurrentNoofKValues = 1;
        end
        
        %define bounds
        lb = zeros(1,CurrentNoofKValues); %lower bound.
        ub = lb+1; %upper bound, may need to be changed for faster reactions
        
        
        %Define the function to optimise:
        f = @(input)Reaction_2(input,KV, Results2, KeepV, Conc,a,CurrentCombinationsFunctions,I,TotalCombinations,Time,QuantifyWhich);
        
        try
            %Completes the timed optimisation of the current reaction model.
            tic
            [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
            toc
            
            %Compile all information into the relevant tables.
            ErrorTable(CombNo,I) = fval;
            KTable{CombNo,I} = num2str(x);
            NoOfKTable(CombNo,I) = length(x);
        catch %in case of ODE integration failure.
            ErrorTable(CombNo,I) = 100;
            KTable{CombNo,I} = '0';
            NoOfKTable(CombNo,I) = 100;
        end
    end
    %Update user on progress:
    sprintf('Optimisation block %d of %d completed.',CombNo,NumberofFunctionArrays)
    %Save current data in case of system failure:
    save('datafile.mat','ErrorTable','KTable','NoOfKTable');
end

%%
%Statistical analysis, all kinetic fitting completed
keyboard

% %Reshape all of the data in order to use AIC.
ErrorTable = reshape(ErrorTable',[],1);
KTable = reshape(KTable',[],1);
NoOfKTable = reshape(NoOfKTable',[],1);

NData = size(Results2,1);
for i = 1:size(ErrorTable,1)
    %Calculate AIC for every reaction model.
    AICraw(i) = NData*(log(ErrorTable(i)/NData))+(2*NoOfKTable(i))+(2*NoOfKTable(i)*((NoOfKTable(i)+1)/(NData-NoOfKTable(i)-1)))+(NData*log(2*pi))+NData;
end
%
% %Sort all of the data in descending ranked order:
[AIC_sorted,AIC_order] = sort(AICraw);
for i = 1:size(ErrorTable,1)
    SortedKTable{i} = KTable{AIC_order(i)};
    SortedErrorTable(i) = ErrorTable(AIC_order(i));
    SortedNoOfKTable(i) = NoOfKTable(AIC_order(i));
end

%Last line of code to place a break to keep all of the workspace:
keyboard;

