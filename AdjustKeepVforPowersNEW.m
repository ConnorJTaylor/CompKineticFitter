function [AdjustedKeepVP,AdjustedKeepVR,AdjKeepVPositions] = AdjustKeepVforPowersNEW(KeepVP,KeepVR,catalytic)
%means you cannot have more than two non integer power reactants
%e.g. A^0.5 + B^0.5 -> C
noninteger2plus = 0;

NumberOfProdReactions = 0;
AdjKeepVPositions = [];
AdjustedKeepV = [];
AdjustedKeepVP = [];
AdjustedKeepVR = [];
powers = [0,0.5,1];
for i = 1:size(KeepVR,2)
    
    if any(KeepVR(:,i)==2) == 1 %if there's a second order reaction, don't compile different orders for that reaction
        Combins = 1;
        CurrentKeepVR = KeepVR(:,i);
        CurrentKeepVP = KeepVP(:,i);
    else
        
        
        NumberofReactants = nnz(KeepVR(:,i)); %nnz is "number of nonzero matrix elements"
        
        if NumberofReactants == 1
            Combins = permn([0,1],NumberofReactants)
        else
        Combins = permn(powers,NumberofReactants);
        end
        KeepWhich = ones(size(Combins,1),1); %Decide on which combins to keep
        for ii = 1:size(Combins,1)
            Currentpower = sum(Combins(ii,:));
            if catalytic == 0
                if Currentpower > 2
                    KeepWhich(ii) = 0;
                elseif Currentpower < 1
                    if NumberofReactants == 2
                        KeepWhich(ii) = 0;
                    end
                end
            else
                if Currentpower > 3 %you're allowed overall third order if it's catalytic
                    KeepWhich(ii) = 0;
                elseif Currentpower < 1
                    if NumberofReactants == 2
                        KeepWhich(ii) = 0;
                    end
                end
            end
        end
        if NumberofReactants == 2
            Combins = Combins(logical(KeepWhich),:);
            KeepWhich = ones(size(Combins,1),1);
            for ii = 1:size(Combins,1)
                %can't have two non-integer powers, for example, k[A]^0.5*[B]^1.5
                if noninteger2plus == 0
                    if ~mod(Combins(ii,1),1) == 0 && ~mod(Combins(ii,2),1) ==0
                        KeepWhich(ii) = 0;
                    end
                end
            end
            Combins = Combins(logical(KeepWhich),:); %final cut to get the desired powers combinations
        end
        
        %create all of the new KeepVR columns
        Combins(Combins==0)=9; %replaces zero values with 9, which is the code for zero order
        CurrentKeepVR = KeepVR(:,i);
        TempCurrentKeepVR = zeros(size(CurrentKeepVR,1),1);
        CurrentKeepVP = KeepVP(:,i);
        ReactantPos = find(CurrentKeepVR~=0);
        CurrentKeepVR = [];
        for ii = 1:size(Combins,1)
            for iii = 1:size(ReactantPos,1)
                TempCurrentKeepVR(ReactantPos(iii)) = Combins(ii,iii);
            end
            CurrentKeepVR = [CurrentKeepVR TempCurrentKeepVR];
        end
    end %end of second order checker
    
    if catalytic(1) == 1
        %also catalytic
        %[catalytic system?, catalytic in which species?]
        
        catorders = [0.5, 1, 2]; %catalytic orders possible
        catalyticspecies = catalytic(2);
        CatalyticCurrentKeepVP = [];
        CatalyticCurrentKeepVR = [];
        for catalystordercount = 1:size(catorders,2) %define the reactants to be modified
            CatalyticCurrentKeepVR = [CatalyticCurrentKeepVR, CurrentKeepVR];
        end
        
        for catalystordercount = 1:size(CatalyticCurrentKeepVR,2) %define the products to be modified
            CatalyticCurrentKeepVP = [CatalyticCurrentKeepVP, KeepVP(:,i)];
        end
        
        
        counter = 1;
        for catalystordercount = 1:size(catorders,2)
            for catcounter= 1:size(CurrentKeepVR,2)
                CatalyticCurrentKeepVP(catalyticspecies,counter) = CatalyticCurrentKeepVP(catalyticspecies,counter)+catorders(catalystordercount);
                if CatalyticCurrentKeepVR(catalyticspecies,counter) == 9
                    CatalyticCurrentKeepVR(catalyticspecies,counter) = 0;
                end
                CatalyticCurrentKeepVR(catalyticspecies,counter) = CatalyticCurrentKeepVR(catalyticspecies,counter)+catorders(catalystordercount);
                
                counter = counter+1;
            end
        end
        
        AdjustedKeepVP = [AdjustedKeepVP CatalyticCurrentKeepVP];
        AdjustedKeepVR = [AdjustedKeepVR CatalyticCurrentKeepVR];
        
        
        TotalnumberKeepV = size(AdjustedKeepVR,2);
        KeepReactionsLogical = ones(1,TotalnumberKeepV);
        TempAdjustedKeepVR = AdjustedKeepVR; %need this to get rid of all 9 orders (zero orders)
        TempAdjustedKeepVR(TempAdjustedKeepVR>8.5)=0;
        
        for i = 1:TotalnumberKeepV
            AnyaboveR = max(TempAdjustedKeepVR(:,i)>2); %any above 2 in reactants?
            AnyaboveP = max(AdjustedKeepVP(:,i)>2); %any above 2 in products?
            
            if AnyaboveR == 1 || AnyaboveP == 1
                KeepReactionsLogical(i) = 0;
            end
            
        end
        
        AdjustedKeepVP = AdjustedKeepVP(:,logical(KeepReactionsLogical));
        AdjustedKeepVR = AdjustedKeepVR(:,logical(KeepReactionsLogical));
        
        
        
        if isempty(AdjKeepVPositions) == 1
            AdjKeepVPositions = [1,size(AdjustedKeepVR,2)];
        else
            lastposition = AdjKeepVPositions(end,end);
            AdjKeepVPositions = [AdjKeepVPositions; lastposition+1 size(AdjustedKeepVR,2)];
        end
        
        
    else 
        
        %remove any reactions where the order is >2 in anything or >2 in product
        TotalnumberKeepV = size(AdjustedKeepVR,2);
        KeepReactionsLogical = ones(1,TotalnumberKeepV);
        TempAdjustedKeepVR = AdjustedKeepVR; %need this to get rid of all 9 orders (zero orders)
        TempAdjustedKeepVR(TempAdjustedKeepVR>8.5)=0;
        
        for i = 1:TotalnumberKeepV
            AnyaboveR = max(TempAdjustedKeepVR(:,i)>2); %any above 2 in reactants?
            AnyaboveP = max(AdjustedKeepVP(:,i)>2); %any above 2 in products?
            
            if AnyaboveR == 1 || AnyaboveP == 1
                KeepReactionsLogical(i) = 0;
            end
            
        end
        
        AdjustedKeepVP = AdjustedKeepVP(:,logical(KeepReactionsLogical));
        AdjustedKeepVR = AdjustedKeepVR(:,logical(KeepReactionsLogical));

        for iv = 1:size(CurrentKeepVR,2)
            AdjustedKeepVP = [AdjustedKeepVP CurrentKeepVP];
        end
        AdjustedKeepVR = [AdjustedKeepVR CurrentKeepVR];
        if isempty(AdjKeepVPositions) == 1
            AdjKeepVPositions = [1,size(CurrentKeepVR,2)];
        else
            lastposition = AdjKeepVPositions(end,end);
            AdjKeepVPositions = [AdjKeepVPositions; lastposition+1 lastposition+size(CurrentKeepVR,2)];
        end
    end
end

end %function end





