function E = TrainingError(W, iter, max_iter)
% Calculates the training error of the Genetic Expression Programming (GEP) model.
% Parameters:
%   - W: The weight vector containing genetic parameters and constants.
%   - iter: Current iteration number.
%   - max_iter: Maximum number of iterations.
% Returns:
%   - E: The training error of the GEP model.

% Global variables
global gep;              % GEP structure
global best_gep;         % Best GEP model
global Emin;             % Minimum error
global BestModel;        % Best model string
global Wdebug;           % Debugging weights
global W_min;            % Minimum weights
global inputs;           % Input data
global outputs;          % Output data
global MSEforGEPDE;      % Mean Squared Error
global head_length;      % Length of the head
global gep_tree_length;  % Length of the GEP tree
global best_constants;   % Best constants
global F;                % Function set
global bestGEPframe;     % Best GEP frame
global x_var;            % Input variable
global c_var;            % Constant variable
global isPrint;          % Global flag indicating whether to print debugging information or not.

Wdebug = W;
% Round the weights to the nearest integer
for i = 1:length(W(1:gep_tree_length))
    W(i) = ceil(W(i));
    Wdebug(i) = W(i);
end
max_arity = 2; % Maximum arity for function nodes
% Extract constants from weight vector
constants = W((gep_tree_length + 1):length(W));
% Create GEP model
gep = gep_icreate(head_length, F, max_arity, x_var, c_var, W);
% Display GEP population if specified
if strcmp(isPrint, 'Yes')
    gep_pprint(gep);
end
% Evaluate GEP model
gep = gep_eval(gep, inputs, constants, outputs);
% Display evaluation information if specified
if strcmp(isPrint, 'Yes')
    fprintf("Iteration          = %d\n", iter);
    fprintf("Simulation Progress= %s\n", strcat("%", string((iter / max_iter) * 100)));
    fprintf("GEP Model          = %s\n", gep.gene.etree.string);
    fprintf("MSE                = %f\n", gep.gene.mse);
end
E = gep.gene.mse; % Training error
% Update best model and constants if error is reduced
if E <= Emin
    BestModel = gep.gene.etree.string;
    Emin = E;
    W_min = W;
    best_constants = W((gep_tree_length + 1):length(W));
    MSEforGEPDE = horzcat(MSEforGEPDE, E);
    bestGEPframe = Wdebug;
    best_gep = gep;
end
% Calculate fitness value
fitness = (1 / (Emin + 1));
% Display best GEP model information if specified
if strcmp(isPrint, 'Yes')
    fprintf("BEST GEP Model     = %s\n", BestModel);
    fprintf("BEST MSE           = %f\n", Emin);
    fprintf("FITNESS            = %f\n", fitness);
    fprintf("BEST CONSTANTS     = {");
    for i = 1:length(best_constants)
        fprintf("c%d=%f; ", i, best_constants(i));
    end
    fprintf("}\n");
    fprintf("%s\n", gep_line_create((4 * (gep.head_length + gep.tail_length + 1) + 30)));
end
end