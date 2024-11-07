function [gep] = gep_eval(gep, inputs, constants, outputs)
% Evaluates the GEP population on the given inputs and constants, computing the mean squared error (MSE) for each individual.
% Parameters:
%   - gep: The Genetic Expression Programming (GEP) population structure.
%   - inputs: Input data matrix.
%   - constants: Array of constants used in the expressions.
%   - outputs: Output data vector.
% Returns:
%   - gep: Updated GEP population structure with MSE evaluated for each individual.

global Emin; % Global variable to store the minimum MSE encountered
[numRows, numCols] = size(inputs); % Get the dimensions of the input matrix
numCons = length(constants); % Get the number of constants
x = {}; % Initializes cell array to store input variables
for i = 1:numCols
    input = inputs(:, i);
    x{i} = input;
    eval(['x' num2str(i) ' = input;']); % Assigns input values to variables x1, x2, ...
end
c = {}; % Initializes cell array to store constants
for i = 1:numCons
    cons = constants(i);
    c{i} = cons;
    eval(['c' num2str(i) ' = cons;']); % Assigns constant values to variables c1, c2, ...
end
gep_model_data = 0; % Initializes the model data
MSE = zeros(1, 1); % Initializes an array to store MSE for each individual
% Evaluates each individual in the population
% for i = 1:gep.psize
    % Evaluating the expression can produce results in complex number format, 
    % which may lead to incorrect calculations.
    gep_model = string(strrep(gep.gene.etree.string, '*', '.*')); % Modify the expression string for evaluation
    gep_model = string(strrep(gep_model, '/', './')); % Modify the expression string for evaluation
    gep_model_data = eval(gep_model); % Evaluate the expression
    % The version using 'subs' is slower.
    % y(i,j)=subs(str2sym(string(gep.gene{i}.etree.string)));
    % disp(gep_model_data);
    if(sum(gep_model_data)==0 | isinf(gep_model_data)| isnan(gep_model_data))
        MSE = 1e+6; % Set a large MSE if the result is invalid
    else
        MSE = mean((outputs - gep_model_data).^2, 'all'); % Compute the mean squared error
    end
    % Check for invalid MSE values
    if MSE < 0 || isinf(MSE) || isnan(MSE)
        MSE = 1e+6; % Set a large MSE if the value is invalid
    end
    % Stores the computed MSE for the individual
    gep.gene.mse = MSE;
% end

% Updates the best MSE and model in the population
% for i = 1:gep.psize
    if gep.gene.mse < Emin
        gep.best.mse = gep.gene.mse;
        gep.best.model = gep.gene.etree.string;
        return; % Returns if a better model is found
    end
% end
end

