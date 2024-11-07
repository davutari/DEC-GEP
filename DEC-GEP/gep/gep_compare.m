function [MSE, MAE, RAE, gep_model_data] = gep_compare(gep_model, inputs, outputs, constants)
% This function compares the performance of a GEP model with given inputs and outputs.
% It evaluates the Mean Squared Error (MSE), Mean Absolute Error (MAE), and Relative Absolute Error (RAE).
% Parameters:
%   - gep_model: The GEP model expression in string format.
%   - inputs: The input data matrix.
%   - outputs: The corresponding output data vector.
%   - constants: A vector containing constant values used in the model.
% Returns:
%   - MSE: The Mean Squared Error between model predictions and actual outputs.
%   - MAE: The Mean Absolute Error between model predictions and actual outputs.
%   - RAE: The Relative Absolute Error between model predictions and actual outputs.
%   - gep_model_data: The model predictions for the given inputs.

[numRows, numCols] = size(inputs);
numCons = length(constants);

% Define symbols for input variables (x) and constants (c).
x = {};
for i = 1:numCols
    input = inputs(:, i);
    x{i} = input;
    eval(['x' num2str(i) ' = input;']);
end

c = {};
for i = 1:numCons
    cons = constants(i);
    c{i} = cons;
    eval(['c' num2str(i) ' = cons;']);
end

% Replace '*' with '.*' and '/' with './' to make the expression compatible with MATLAB evaluation.
gep_model = string(strrep(gep_model, '*', '.*'));
gep_model = string(strrep(gep_model, '/', './'));

% Evaluate the GEP model expression using the input data.
gep_model_data = eval(gep_model);

% Calculate evaluation metrics.
MAE = mean(abs(outputs - (gep_model_data)), 'all'); % Mean Absolute Error
MSE = mean((outputs - gep_model_data).^2, 'all'); % Mean Squared Error
RAE = mean(abs(outputs - gep_model_data) ./ abs(outputs), 'all'); % Relative Absolute Error
end