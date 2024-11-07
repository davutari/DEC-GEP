% DEC-GEP (Differential Evolution Chromosomal Gene Expression Programming) 
% tool This tool was developed for data-driven symbolic regression application.   
% DEC-GEP benefits from the advantages of optimizing a complete genotype of 
% expression trees for optimal model extraction from a training datasets 
% This tool implements Differential Evolution (DE) optimization 
% in order to search for optimal expression tree that expresses 
% optimal symbolic regression model from dataset. 

% For detailed information, 
% you can refer to the https://www.sciencedirect.com/science/article/pii/S1568494623001114.

% Citation: Ari, D., & Alagoz, B. B. (2023). A differential evolutionary chromosomal gene expression 
% programming technique for electronic nose applications. Applied Soft Computing, 136, 110093.
% https://doi.org/10.1016/j.asoc.2023.110093

clear all;close all;
clc;
% Global variable declaration of the DEC-GEP tool
% Please do not change global variable declarations
global inputs;            % Input data for symbolic regression model
global outputs;           % Output data for symbolic regression model
global head_length;       % Length of Head section in the GEP chromosome
global best_constants;    % The optimal set of constants after the DEC-GEP algorithm is completed.
global BestModel;         % The optimal symbolic model after the DEC-GEP algorithm is completed.
global F;                 % Function set (F) to defines functions and operators
global isPrint;           % Flag to enable to print report of DEC_GEP at each iteration
global x_var;             % Variable set to define input terminals
global c_var;             % Constant set to define constant terminals
% Sets isPrint flag to report significant parameters at each iteration ('Yes' or 'No')
isPrint = 'No';
% Loads data from the specified file
load('./data/AirQualitySensor');
% Assigns input data from the loaded dataset
inputs = x';
% Assigns output data from the loaded dataset
outputs = yCO';
% Sets the head length of the GEP chromosome
head_length = 16;
% Defines the lower bound for constant values
ConsValLow = -1000;
% Defines the upper bound for constant values
ConsValMax = 1000;


% Defines the function set (F) for the GEP algorithm
% Format of set is the {'Function',number of arguments}
% Each row represents a function and its arity (number of arguments (operands))
% Example:
% F = {'*', 2;   % Multiplication operator with arity 2
%      '+', 2;   % Addition operator with arity 2
%      '-', 2;   % Subtraction operator with arity 2
%      'sin', 1}; % Sinus function with the arity 1

F={'div',2;'ReelSqrt',1;'ReelLog',1;'cos',1;'tan',1;'cot',1;'sin',1;'*',2;'+',2;'-',2};

% Defines the variable terminal set (V)
% This set includes the input variables for the GEP model
x_var = {'x1','x2','x3','x4','x5','x6','x7','x8'};  % Terminal variables vector
% Defines the constant terminal set (C)
% This set includes the constant values that can be used in the GEP model  
c_var = {'c1','c2','c3','c4','c5','c6','c7','c8'};

%In Genetic Programming (GP), the term "arity" refers to the number of arguments 
% or operands that a function or operator takes.
% max_arity: the maximum arity (number of arguments) of the functions in the function set
max_arity=2;

%-----------Differential Evolution(DE) Parameters-----------------------------------
% Maximum Number of Iterations for the DE algorithm
Max_iteration = 5000;
% Population Size for the DE algorithm
SearchAgents_no = 40; 
% Lower Bound of the Scaling Factor in DE
beta_min = 0.2;   
% Upper Bound of the Scaling Factor in DE
beta_max = 1.5;   
% Crossover Probability for the DE algorithm
pCR = 0.2;       
% Message flag for optimization process
% 0: stop messages
% 1: enables messages from the optimization process
MesFlag = 1;
% Enables to plot convergence characteristic (Cost-Iteration curve)
% 0: Does not draw Cost-Iteration plot
% 1: Draws Cost-Iteration plot
IsDraw = 1; 

% Run the DEC_DE algorithm 
DEC_GEP(ConsValLow,ConsValMax,max_arity,Max_iteration,SearchAgents_no, beta_min,beta_max,pCR,MesFlag,IsDraw)

% Prints a separator line for better readability in the command window
fprintf('********************************************************************\n');
% Calculates the predicted data and evaluate performance of the symbolic regression model
[MSE, MAE, RAE, gep_model_data] = gep_compare(BestModel, inputs, outputs, best_constants);
% Prints the Mean Absolute Error (MAE)
fprintf('MAE: %d \n', MAE);
% Prints the Mean Squared Error (MSE)
fprintf('MSE: %d \n', MSE);
% Prints the Relative Absolute Error (RAE)
fprintf('RAE: %d \n', RAE);
% Prints the symbolic regression model
fprintf('Symbolic Model:');
disp(BestModel);
% Checks if the constants should be printed
if (strcmp(isPrint, 'No')) 
    % Print each constant
    for i = 1:length(best_constants)
        fprintf("c%d=%f; ", i, best_constants(i));
    end
end
% Prints a newline for better readability
fprintf("\n");
% Opens figure (2) for plotting the real data vs. predicted data
figure(2);
% Plots the real data and the GEP prediction
plot(1:length(outputs), outputs, 1:length(outputs), gep_model_data, 'LineWidth', 1.75);
xlabel('Data Instance');
ylabel('Value');
% Adds a legend to the plot
legend('Real Data', 'DEC-GEP Estimation');
grid
% Get the current date and time in the specified format and time zone
currentDate = datestr(datetime('now', 'TimeZone', 'Europe/Istanbul'), 'dd-mm-yyyy');
currentTime = datestr(datetime('now', 'TimeZone', 'Europe/Istanbul'), 'HH:MM:SS');
% Replaces colons with dashes in the time string
currentTime = strrep(currentTime, ':', '-');
% Combines the date and time into a single string
SDate = strcat(currentDate, '_', currentTime);
% Creates a filename for saving the workspace
strName = strcat("Workspace_", SDate);
% Saves the current workspace to a file with a name that includes date and
% time
save(strName);