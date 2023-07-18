%Author: Godwin Mruma Gadiel
%Description: Tabu search based algorithm for reproducing the results for
%the manuscript tittled "Energy Efficient Phase Interpolator based Hybrid
%Beamforming Architecture for massive MIMO System.
function [FRF,WRF, Rate] = functionTS_IP_ET_TX_TX(H,P,Ns,At,Ar,NRF,SNR,Fopt)
% Tabu search algorithm
[~,Nt] = size(H);
[~,N] = size(At);

%%% Problem definition
CostFunction = @ function_information_rate;  % Cost Function

%Number of actions is defined by number of antennas
% Note that the total number of action is Nt*U however U times will be
% applied inside the loop
nAction1 = Nt;

%%% Tabu Search Parameters

MaxIt = 60;                      % Maximum Number of Iterations

%TL1 = 5*nAction1;   % in our case we considered the tabu list to be infinity
 
%%% Initialization

% Create Empty Individual Structure
empty_individual.analogPrecoder = [];
empty_individual.analogCombiner = [];
empty_individual.Systemrate = [];
empty_individual.SolIndex_tx = [];
empty_individual.SolIndex_rx = [];

% Create Initial Solution
sol = empty_individual;
%I_solgen_tx = function_solgen_in(N,Ns);  % generating initial solution for tx
[ F_RF, F_BB, Sol_ind_tx ] = functionOMP( Fopt,Ns, NRF, At );
I_solgen_tx = zeros(1,N);
I_solgen_tx(1,Sol_ind_tx) = 1;
sol.SolIndex_tx = Sol_ind_tx;
%I_solgen_rx = function_solgen_in(N,Ns);  % generating initial solution for rx
[~, ~, Sol_ind_rx] = functionOMPSPMMSE(H,F_RF, F_BB,Ar,Ns,SNR);
I_solgen_rx = zeros(1,N);
I_solgen_rx(1,Sol_ind_rx) = 1;
sol.SolIndex_rx = Sol_ind_rx; 
sol.ananogPrecoder = At(:,sol.SolIndex_tx);
sol.analogCombiner = Ar(:,sol.SolIndex_rx);
sol.Systemrate = CostFunction(H, P, Ns, sol.ananogPrecoder, sol.analogCombiner);

% Initialize Best Solution Ever Found
BestSol = sol;

% Array to Hold Best Costs
BestCost = zeros(MaxIt,1);

% Initialize Action Tabu Counters
TC_tx = {};
TC_rx = {};


%%% Tabu Search Main Loop
cn = 0;   % initialize counter for early termination
for it =1:MaxIt
    
    bestnewsol.Systemrate = -inf;
    
    % Apply Actions
    for i = 1:nAction1
        % this action is for the first iteration(just for initialization)
        if it == 1
            bestnewsol.SolIndex_tx = sol.SolIndex_tx;
            bestnewsol.SolIndex_rx = sol.SolIndex_rx;
            %bestnewsol.SolIndex_rx = I_solgen_rx;
        end
        I_solgen_new_tx = function_solgen(N,NRF,bestnewsol.SolIndex_tx);  % generating solution for tx
        I_solgen_new_rx = function_solgen(N,Ns,bestnewsol.SolIndex_rx);  % generating solution for rx
        LS_tx = zeros(1,length(TC_tx));
        LS_rx = zeros(1,length(TC_rx));
        for k = 1:length(TC_tx)

            LS_tx(k) = sum(ismember(I_solgen_new_tx,TC_tx{k}));
            LS_rx(k) = sum(ismember(I_solgen_new_tx,TC_rx{k}));
        end
                if ismember(Ns,LS_tx) ~= 1
                    newsol.SolIndex_tx = I_solgen_new_tx;
                    newsol.SolIndex_rx = I_solgen_new_rx;
                    newsol.ananogPrecoder = At(:,newsol.SolIndex_tx);
                    newsol.analogCombiner = Ar(:, newsol.SolIndex_rx);
                    newsol.Systemrate = CostFunction(H, P, Ns, newsol.ananogPrecoder, newsol.analogCombiner);
                   


                    newsol.ActionIndex = {i};
                else
                    newsol = sol;             % keep the current solution
                    newsol.ActionIndex = {i};
                end

                if newsol.Systemrate  > bestnewsol.Systemrate 
                    bestnewsol = newsol;
                end
            
        
    end
   
    % Update Current Solution
    sol = bestnewsol;
    
    % Update Tabu List*** We assume tabu list is too large and therefore
   TC_tx{i} = bestnewsol.SolIndex_tx;
   TC_rx{i} = bestnewsol.SolIndex_rx;
    % Update Best Solution Ever Found
    if sol.Systemrate > BestSol.Systemrate
        BestSol = sol;
    end
    
    % Save Best Cost Ever Found
    BestCost(it) = BestSol.Systemrate;
      
    
    % early termination
     if it > 1
        if BestCost(it-1) == BestCost(it)
            cn = cn + 1;
            if cn == 10
                break
            end
        else
            cn = 0;
        end
        
    end
    
end  % end of iteration
F = BestSol.SolIndex_tx;  % the best selected indices
W = BestSol.SolIndex_rx;
FRF = At(:,F);
WRF = Ar(:,W);
Rate = BestSol.Systemrate;


end% end of the function

