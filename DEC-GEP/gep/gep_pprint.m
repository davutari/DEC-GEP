function gep_pprint(gep)
% Pretty-print the genotype of each individual in the GEP population.
% Parameters:
%   - gep: The Genetic Expression Programming (GEP) population structure.
% Displays:
%   - The genotype of each individual in the population.
% Prints a line separator
fprintf("%s\n", gep_line_create((4 * (gep.head_length + gep.tail_length + 1) + 30)));
    fprintf("Genotype {"); % Start of genotype representation
    % Prints each gene in the individual's genotype
    for j = 1:length(gep.gene.etree.gene)
        fprintf(" %s", string(gep.gene.etree.gene{j}));
    end
    fprintf("}\n"); % End of genotype representation
end