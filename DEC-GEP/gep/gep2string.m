function [stree] = gep2string(etree)
% Convert a genotype tree structure into a string representation using Karva Notation.
% Parameters:
%   - etree: The genotype tree structure to be converted.
% Returns:
%   - stree: The string representation of the genotype tree.
%genotip{+  x2 *  *  *  -  Q  Q  c1 x2 x2 x1 c1 c2 x2 x2 x1 x1 x1 x1 c2 c2 c2 x2 c2 x1 c2 c2 x2 c2 c1 c1 c2 }
%fenotip y=(x2+(((x2-x2)*Q(x1))*(Q(c1)*c1)))
%The genotype in string array format is converted into matrix format according to Karva Notation (NX3)
% --ETree Index Matrix Form---
% |Node|Node Left|Node Right|
% |1   |2        |3         |
% |2   |NaN      |NaN       |
% |3   |4        |5         |
% |4   |6        |7         |
% |5   |8        |9         |
% |6   |10       |11        |
% |7   |12       |NaN       |
% |8   |13       |NaN       |
% |9   |NaN      |NaN       |
% |10  |NaN      |NaN       |
% |11  |NaN      |NaN       |
% |12  |NaN      |NaN       |
% |13  |NaN      |NaN       |
% |14  |NaN      |NaN       |
% |15  |NaN      |NaN       |
% |16  |NaN      |NaN       |
% |17  |NaN      |NaN       |
% |18  |NaN      |NaN       |
% |19  |NaN      |NaN       |
global usedlist;
% Initialize a list to keep track of used indices.
usedlist = zeros(1, length(etree.gene));
% Initialize a matrix to represent the genotype tree structure.
matrix_form = zeros(length(etree.gene), 3);
% Populate the matrix with indices and connections between nodes.
for i = 1:length(etree.gene)
    matrix_form(i, 1) = i;
    if (etree.itemtype(i) == 1) % If the item is a function.
        matrix_form(i, 2) = getusedindex(); % Assign an index for the left node.
        if (etree.itemport(i) == 2)
            matrix_form(i, 3) = getusedindex(); % Assign an index for the right node.
        else
            matrix_form(i, 3) = NaN; % If the function is unary, mark the right node as NaN.
        end
    else % If the item is a terminal.
        matrix_form(i, 2) = NaN; % Mark both left and right nodes as NaN.
        matrix_form(i, 3) = NaN;
    end
end
% Convert the genotype tree structure into a string using the matrix2string function.
stree = matrix2string(etree, 1, matrix_form);
end
function [sout] = matrix2string(etree, index, matrix_form)
    % This function converts the genotype tree (etree) into a string expression.
    % It traverses the tree recursively starting from the given index.
    % Parameters:
    %   - etree: The genotype tree structure containing genes and their connections.
    %   - index: The index of the current node in the genotype tree.
    %   - matrix_form: The matrix representation of the genotype tree.
    % Returns:
    %   - sout: The string representation of the subtree rooted at the given index.
    if (index <= length(etree.gene))
        % Check if the index is within the bounds of the gene array.
        if (etree.itemtype(index) == 1)
            % If the current node is a function node:
            if (etree.itemport(index) == 1)
                % If it's a unary operation (arity 1):
                % Recursively call matrix2string for the left subtree.
                t = matrix2string(etree, matrix_form(index, 2), matrix_form);
                % Construct the string expression with the current function and its argument.
                sout = strcat(etree.gene{index}, '(', t, ')');
            else
                % If it's a binary operation (arity 2):
                % Recursively call matrix2string for the left and right subtrees.
                sleft = matrix2string(etree, matrix_form(index, 2), matrix_form);
                sright = matrix2string(etree, matrix_form(index, 3), matrix_form);
                if (isselffunction(etree.gene{index}) == 1)
                    % If the current function is self-referencing (e.g., +, *, etc.),
                    % Wrap the left and right subtree expressions with parentheses.
                    sout = strcat(etree.gene{index}, '(', '(', sleft, ')', ',', '(', sright, ')', ')');
                else
                    % Otherwise, construct the string expression normally.
                    sout = strcat('(', sleft, ')', etree.gene{index}, '(', sright, ')');
                end
            end
        else
            % If the current node is a terminal node (a variable or constant),
            % Simply return its corresponding gene from the gene array.
            sout = etree.gene{index};                 
        end
    else
        % If the index is out of bounds, return an empty string.
        sout = etree.gene{index};
    end
end
function [uindex] = getusedindex()
    % Finds an unused index in the 'usedlist' array.
    % Parameters:
    %   - None
    % Global Variables:
    %   - usedlist: Array indicating the usage status of indices.
    % Output:
    %   - uindex: The index in 'usedlist' where the value is 0. Returns -1 if no unused index is found.
    global usedlist;
        for i = 2:length(usedlist)
            if (usedlist(1,i) == 0)
                usedlist(1,i) = 1;
                uindex = i;
                return
            end
        end
        uindex = -1;
end

function [isselffunc] = isselffunction(func)
    % Checks if the given function is a basic MATLAB operator.
    % Parameters:
    %   - func: The function name to be checked.
    % Returns:
    %   - isselffunc: Returns 1 if the function is a basic MATLAB operator, otherwise returns 0.
    % Global Variable:
    %   - matlab_op: String containing basic MATLAB operators: '+', '/', '*', '-'.
    matlab_op='+, /, *, -';%Matlab operatörleri
    isselffunc = contains(matlab_op,string(func));
    if(isselffunc==0)
        isselffunc =1;
    else
        isselffunc =0;
    end
end