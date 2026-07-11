%% Copyrigth Sidhant Manoj Barai


clear all;
clc;

%% Sampling Plan
% Initialising gains
initial_gains = [10, 10]; % Kp = 10 and KI = 10; gains

% 1st Sampling - FullFactorial (X - Sampling plan for Z)
ffsampling = fullfactorial(initial_gains, 1); 
% Default = 1 for equal spacing of samples near the edges
ffsampling = 10*ffsampling + eps; % Scaling

% 2nd Sampling - Rlh Sampling (X - Sampling plan for Z)
rlhsampling = rlh(100, 2, 3);
rlhsampling = 10*rlhsampling + eps; % Scaling

% 3rd Sampling Methods - Sobol (X - Sampling plan for Z)
P = sobolset(2);
sobolsampling = net(P, 100);
sobolsampling = 10*sobolsampling + eps; % Scaling

%% Scatter plots for the sampling of KP and KI gains
figure;
scatter(ffsampling(:, 1), ffsampling(:, 2), 'magenta')
xlabel('K_p');
ylabel('K_i');
title('Full Factorial Sampling')

figure; 
scatter(rlhsampling(:, 1), rlhsampling(:, 2), 'magenta')
xlabel('K_p');
ylabel('K_i');
title('Random Latin Hypercube Sampling')

figure;
scatter(sobolsampling(:, 1), sobolsampling(:, 2), 'magenta')
xlabel('K_p');
ylabel('K_i');
title('Sobol Sampling')

%% Assessing Sampling Plans using Metrics

ff_Metric = mmphi(ffsampling, 1, 2); % Full Factorial Sampling Metric
rlh_Metric = mmphi(rlhsampling, 1, 2); % Rlh Sampling Metric
sobol_Metric = mmphi(sobolsampling, 1, 2); % Sobol Metric

% Selection of best metric and using the sampling method for future use
% Displaying the name of the best sampling method
if ff_Metric < rlh_Metric && ff_Metric < sobol_Metric
    bestSampling = ffsampling;
    disp('Full Factorial Sampling is the best sampling method')
elseif rlh_Metric < ff_Metric && rlh_Metric < sobol_Metric
    bestSampling = rlhsampling;
    disp('Rlh sampling is the best sampling method')
else
    bestSampling = sobolsampling;
    disp('Sobol Sampling is the best sampling method')
end

% Evaluate Design to compare relationship between 
% performance criteria and design
Z = evaluateControlSystem(bestSampling);

% Plotting the correlation 
figure;
label = {'Largest CL Pole', 'Gain Margin', 'Phase Margin', 'RiseTime', ...
    'Peak Time', 'MaxOShoot', 'MaxUShoot', 'SettlingTime', 'SS Error',...
    'ControlEffort'};
[h, ax] = plotmatrix(Z);

% Setting the labels
numVars = size(Z, 2);
for i = 1:numVars
    ax(numVars, i).XLabel.String = label{i};
    ax(i, 1).YLabel.String = label{i};
end

% Plotting to understand the relationship between performance criteria 
% and design using parallel plotting
figure; 
parallelplot(array2table(Z, 'VariableNames', label));
title('Parallel Plot of Performance Criteria');

%% Optimization
% Changing Gain Margin and Phase Margin as postprocessing. This is done 
% before the actual optimization process
Z_value = optimizeControlSystem(bestSampling);

% %Initializaion for hypervolume indicator
% HV = [];
% reference = [max(Z_value)];

% Optimization using NSGA-II technique for 250 iterations
i = 1;
sample_plan = bestSampling;
while i<51
    
    disp(i);

    % Saving Parent plan
    parent_sampling_plan = sample_plan;
    
    % Pre-Defining the Priority Levels:
    % Hard Constraint = 3; High = 2; Moderate = 1; Low = 0;
    if i == 1 && i<=50
        priority = [3 2 2 1 0 1 0 0 1 2];
        goals = [1 -6 20 2 10 10 8 20 1 0.67];
    end

    % Incorporating rank preferability
    [rank, ClassV] = rank_prf(Z_value, goals, priority);

    % Crowd Sorting Method
    crowd = crowding(Z_value, rank);

    % Inverting rankings to make sure the btwr function considers 
    % fitness = 3 as the better rank rather than fitness = 0
    rank = max(rank) - rank;

    % Binary Tournament Selection using the fitness values
    selection = btwr([rank, crowd], 100);

    % Defining boundary conditions for mutation and cross-over
    boundary = [0 0; 1 1]; % Assumptions

    % Cross-Over and Mutation variation
    crossOver = sbx(parent_sampling_plan(selection,:), boundary);

    % Mutation to get the child after variation
    child = polymut(crossOver, boundary);

    % Joining the child and parent
    parent_child = [parent_sampling_plan; child];

    % Evaluating the parent-child concatenation 
    Z_value_PC = optimizeControlSystem(parent_child);

    % Doing the same ranking preferability for the parent-child combo
    [rank_PC, ClassV] = rank_prf(Z_value_PC, goals, priority);

    % Crowd Sorting for the parent-child combo
    crowd_PC = crowding(Z_value_PC, rank_PC);

    % Using NSGA-II selection for next iteration
    nsga_II = reducerNSGA_II(parent_child, rank_PC, crowd_PC, 100);

    % Finding the new NSGA-II reduced population
    parent_child_new = parent_child(nsga_II,: );

    % Evaluating the new reduced population
    Z_value = optimizeControlSystem(parent_child_new);

    % % HyperVolume Indicator
    % HV_indicator = log10(Hypervolume_MEX(Z_value, reference));
    % % Updating HyperVolume Vector
    % HV = [HV; HV_indicator];

    % Updating the new Sampling Plan
    sample_plan = parent_child_new;

    i = i+1;
end


%% Plotting optimized result of design
figure;
scatter(parent_sampling_plan(:, 1), parent_sampling_plan(:, 2))
xlabel('K_p of new Parent');
ylabel('K_i of new Parent');
title('Best Sampling')

% Again plotting the correlation but this time using the optimized Z value
figure;
[h, ax] = plotmatrix(Z_value);

% Setting the labels
numVars = size(Z_value, 2);
for i = 1:numVars
    ax(numVars, i).XLabel.String = label{i};
    ax(i, 1).YLabel.String = label{i};
end

% Plotting to understand the relationship between performance criteria and
% design using parallel plotting after optimization
figure; 
parallelplot(array2table(Z_value, 'VariableNames', label));
title('Parallel Plot of Performance Criteria');

%% Plotting HyperVolume
figure;    
plot(HV, 'LineWidth', 2, 'Color', 'r');
    xlim([0, 51]);
    ylim([208, 209]);
    title('Hypervolume Convergence plot');
    xlabel('Iteration');
    ylabel('Hypervolume (log10 scale)');

%% Function for setting up gain margin and phase margin
% evaluation for more samples

function Z = optimizeControlSystem(P)
    Z = evaluateControlSystem(P);
    % Gain Margin change
    Z(:,2) = -20*log10(Z(:,2));
    % Phase Margin change
    Z(:, 3) = abs(Z(:, 3) - 50);
    % Remove inf values
    Z(isinf(Z)) = 1e6;
end
