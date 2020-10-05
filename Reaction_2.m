function [Error] = Reaction_2(input,KV, Results2, KeepV, Conc,a,CurrentCombinationsFunctions,I,TotalCombinations,Time,QuantifyWhich)
%Function to be optimised that calls ODE solver and outputs 'Error' as the
%variable to be optimised. K values are iterated to find the minimum sum of
%squared error between the simulated data and the experimental data.

KV = input; %assign inputs to k values

    %runs the ODE function
    interupt_time = 100; %Maximum time in seconds for the ODE solver to run.
    
    %Sets the output function to incorporate interuption time:
    outputFun = @(Time,Conc,flag)interuptFunODE(Time,Conc,flag,interupt_time);
    opts1 = odeset('OutputFcn',outputFun);
    opts2 = odeset('OutputFcn',outputFun,'Nonnegative',1);

    warning off
    
    %This is the maximum allowed total molar value of all species, this is
    %defined so that ODEs don't form a 'runaway integration' where values
    %quickly escalate by orders of magnitude.
    TotalConc = sum(Conc)*10;
   
    
    %ODE solver:
    [TimeData,ConcData] = ode15s(@(Time,Conc)ReactionKineticLaws_2(Time,Conc,KV,a,CurrentCombinationsFunctions,I,TotalConc),Time,Conc,opts2);

    %Assign quantified species:
    ConcData2 = [];
    for i = 1:size(QuantifyWhich,2)
        ConcData2 = [ConcData2 ConcData(QuantifyWhich(i),:)];
    end
    
    Error = sse(ConcData2-Results2);
    %Sum of squared error comparison, serves as the best comparison to
    %experimental results after the particular run using the particular K
    %inputs.
end





