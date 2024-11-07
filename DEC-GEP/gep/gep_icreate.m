function gep = gep_icreate(head_length, functions, max_arity, variables, constants, W)
% Creates an initial population of Genetic Expression Programming (GEP) individuals.
% Parameters:
%   - psize: Population size.
%   - head_length: Length of the head.
%   - functions: Set of functions available for the GEP.
%   - max_arity: Maximum arity of functions.
%   - variables: Set of variables available for the GEP.
%   - constants: Set of constants available for the GEP.
%   - W: Weight vector containing genetic parameters.
% Returns:
%   - gep: Initial population of GEP individuals.
% Initializes GEP structure
gep.generation = 1;                          % Current generation
gep.head_length = head_length;               % Length of the head
gep.tail_length = head_length * (max_arity - 1) + 1; % Length of the tail
gep.functions = functions;                   % Set of functions
gep.variables = variables;                   % Set of variables
gep.constants = constants;                   % Set of constants
gep.best.model = NaN;                        % Best model placeholder
gep.best.mse = NaN;                          % Best mean squared error placeholder
gep.gene.fitness = 0;                 % Initialize fitness
gep.gene.mse = 0;                      % Initialize mean squared error
gep.gene.etree = de2gep(gep, W);       % Create expression tree for the individual

end
