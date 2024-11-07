% DEC-GEP (Differential Evolution Chromosomal Gene Expression Programming) improves 
% the efficiency of GEP chromosomes through a structured design. The DEC-GEP chromosome starts with a start codon, 
% ensuring valid expression tree initiation, followed by a head root (function-only codes), 
% head body (functions, variables, and constants), 
% tail (variables and constants), and a numeric box for constant values.
% This structure allows DE to optimize the entire chromosome, enhancing the 
% exploration and exploitation capabilities of GEP. The method integrates numerical optimization within the evolutionary process,
% addressing traditional GEP limitations.
% For detailed information, you can refer to the https://www.sciencedirect.com/science/article/pii/S1568494623001114.

% Citation: Ari, D., & Alagoz, B. B. (2023). A differential evolutionary chromosomal gene expression 
% programming technique for electronic nose applications. Applied Soft Computing, 136, 110093.
% https://doi.org/10.1016/j.asoc.2023.110093

function DEC_GEP(ConsValLow,ConsValMax,max_arity,Max_iteration,SearchAgents_no, beta_min,beta_max,pCR,MesFlag,IsDraw)

% Declare global variables used in the DE-GEP algorithm
global MSEforGEPDE;       % Mean Squared Error for GEP-DE
global Emin;              % Minimum error threshold
global inputs;            % Input data for the model
global outputs;           % Output data for the model
global NumofCons;         % Number of constants used in the model
global head_length;       % Head length of GEP chromosome
global gep_tree_length;   % Length of the GEP tree
global best_constants;    % Best set of constants found
global BestModel;         % Best model found
global F;                 % Function set (F)
global isPrint;           % Flag to determine if constants should be printed
global x_var;             % Variable for input data
global c_var;             % Variable for constants

% Determine the number of rows and columns in the input data
[numRows, numCols] = size(inputs);
% Number of variables (input features) in the dataset
NumofVar = numCols;
% % Set the head length of the GEP chromosome
% head_length = 16;
% % Define the lower bound for constant values
% ConsValLow = -1000;
% % Define the upper bound for constant values
% ConsValMax = 1000;
% Initialize the MSE for GEP-DE to a large value
MSEforGEPDE = 10^10;
% Initialize the minimum error threshold to a large value
Emin = 10^10;
% Add custom functions path 
addpath(string(strcat(pwd,'\functions')));
% Add GEP functions path 
addpath(string(strcat(pwd,'\gep')));
% Adding the path with optimization methods such as DE, PSO, GWO
addpath(string(strcat(pwd,'\optimization methods')));
% Add data path 
addpath(string(strcat(pwd,'\data')));
%Operator list={'Function',port}
% Define the function set (F) for the GEP algorithm
% Each row represents a function and its arity (number of arguments)
% F = {'*', 2;   % Multiplication function with arity 2
%      '+', 2;   % Addition function with arity 2
%      '-', 2;
%      'sin', 1};
     %'exp', 1};  % Subtraction function with arity 2
% Alternative function set including division
% F = {'div', 2; '*', 2; '+', 2; '-', 2};
F={'div',2;'ReelSqrt',1;'ReelLog',1;'cos',1;'tan',1;'cot',1;'sin',1;'*',2;'+',2;'-',2};
% Define the terminal variable set (x_var)
% This set includes the input variables for the GEP model
x_var = {'x1','x2','x3','x4','x5','x6','x7','x8'};  % Terminal variables vector ['x1', 'x2']
% Define the constant variable set (c_var)
% This set includes the constant values that can be used in the GEP model
%c_var = {'x1','x2','x3','x4','x5','x6','x7','x8'};  % Constant variables vector ['c1', 'c2', 'c3', 'c4', 'c5']max_arity = 2;%In Genetic Programming (GP), the term "arity" refers to the number of arguments or operands that a function or operator takes. 
c_var = {'c1','c2','c3','c4','c5','c6','c7','c8'};
%c_var = {};
% Calculate the length of the tail in a Gene Expression Programming chromosome
% The tail length is calculated based on the head length and the maximum arity
% head_length: the length of the head of the chromosome
% max_arity=2;
% max_arity: the maximum arity (number of arguments) of the functions in the function set
% The tail length is computed using the formula:
% tail_length = head_length * (max_arity - 1) + 1
% This ensures that the tail is long enough to accommodate the arguments required
% by the functions in the function set, allowing for valid expressions to be generated.
% Calculate the tail length of the GEP chromosome
% Tail length is determined by the formula: head_length * (max_arity - 1) + 1
tail_length = head_length * (max_arity - 1) + 1;
% Calculate the total length of the GEP expression tree (chromosome)
% This is the sum of the head length and tail length
gep_tree_length = head_length + tail_length;
% Determine the number of constant variables in the GEP model
% NumofCons is the length of the constant variable array (c_var)
NumofCons = length(c_var);
% Calculate the total parameter size for the GEP model
% This includes the size of the function set (F), number of variables (NumofVar), and number of constants (NumofCons)
param_size = length(F) + NumofVar + NumofCons;
%-----------DE (Differential Evolution) parameters-----------------------------------
% % Maximum Number of Iterations for the DE algorithm
% Max_iteration = 5000;
% % Population Size for the DE algorithm
% SearchAgents_no = 40; 
% % Lower Bound of the Scaling Factor for mutation in DE
% beta_min = 0.2;   
% % Upper Bound of the Scaling Factor for mutation in DE
% beta_max = 1.5;   
% % Crossover Probability for the DE algorithm
% pCR = 0.2;       
% % Message flag for optimization process
% % 0: stop messages
% % 1: enables messages from the optimization process
% MesFlag = 1; 
% Configuration for Optimization Algorithm
% Objective function handle to calculate the training error
fobj = @(W, iter, max_iter) TrainingError(W, iter, max_iter);  
% Lower and upper bounds for the root element to ensure it is a function
root_lb = [0]; % Ensure the first element is from the function set
root_ub = [length(F)];
% Lower and upper bounds for the head section (excluding root) of the GEP chromosome
head_lb = ones(1, head_length - 1); % Head section (excluding root) can take any element
head_ub = param_size * ones(1, head_length - 1); 
% Lower and upper bounds for the tail section of the GEP chromosome
tail_lb = (length(F) + 1) * ones(1, tail_length); % Tail section can only take variables and constants
tail_ub = param_size * ones(1, tail_length); 
% Lower and upper bounds for the constant values
ConsValue_lp = ConsValLow * ones(1, NumofCons); % Lower bound for constant values
ConsValue_ub = ConsValMax * ones(1, NumofCons); % Upper bound for constant values
% Combine all lower bounds into one vector
lb = [root_lb head_lb tail_lb ConsValue_lp]; 
% Combine all upper bounds into one vector
ub = [root_ub head_ub tail_ub ConsValue_ub]; 
% Dimension of the optimization problem
dim = length(lb); 
% Execute the DE algorithm to optimize the GEP model
% Outputs are the best solution (BestSol), the best cost (BestCosts), and the cost history (BestCostHist)
[BestSol, BestCosts, BestCostHist] = de(fobj, dim, lb, ub, Max_iteration, SearchAgents_no, beta_min, beta_max, pCR, MesFlag);
%-------------------------------------------------------
% Create a figure for plotting the best cost history
if IsDraw==1
figure(1);
% Plot the history of the best cost over iterations
plot(BestCostHist, 'LineWidth', 1.75);
% Label the x-axis as 'Iteration'
xlabel('Iteration');
% Add a legend to the plot
ylabel('Cost');
grid
end
