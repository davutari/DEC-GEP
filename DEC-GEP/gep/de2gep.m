function [etree] = de2gep(gep, de_array)
  % Ensures the head function has a valid index (at least 1)
  if (de_array(1) <= 0)
    de_array(1) = 1;
  end
  % Extracts the head function from the function set based on its index
  etree.head = gep.functions(de_array(1));
  % Converts the head function name to a string and store it in the gene expression
  etree.gene = strcat(string(etree.head));
  % Sets the item type for the head (function) to 1
  etree.itemtype(1) = 1;
  % Extracts the number of arguments (ports) for the head function
  etree.itemport(1) = cell2mat(gep.functions(de_array(1), 2));
  % Loop through remaining elements in the DE array (excluding head)
  for i = 2:length(de_array(1:(gep.head_length + gep.tail_length)))
    % Gets details (type, name, ports) of the current item based on its index in the DE array
    treeitem = get_items(gep, i, de_array);
    % Stores the item type (function, variable, or constant)
    etree.itemtype(i) = treeitem.type;
    % Appends the item name (function name, variable name, or constant value) to the gene expression
    etree.gene = horzcat(etree.gene, string(treeitem.name));
    % Stores the number of arguments (ports) for the current item
    etree.itemport(i) = treeitem.ports;
  end
  % Converts the expression tree (etree) into a string representation
  etree.string = gep2string(etree);
end
%------------------------------------------------------------------
function [item] = get_items(gep, index, de_array)
    % This function retrieves an item (function, variable, or constant) from the GEP.
    % It looks up the item based on the given index and decision array.
    % Parameters:
    %   - gep: The GEP structure containing functions, variables, and constants.
    %   - index: The index of the item in the decision array.
    %   - de_array: The decision array specifying which items to select.
    % Returns:
    %   - item: A structure containing information about the selected item.
    if de_array(index) <= length(gep.functions)
        % If the index corresponds to a function:
        item.name = gep.functions(de_array(index)); % Get the function name.
        item.type = 1; % Type 1 indicates a function.
        item.ports = cell2mat(gep.functions(de_array(index), 2)); % Get the number of ports.
    elseif ((length(gep.functions) < de_array(index)) && (de_array(index) <= length(gep.functions) + length(gep.variables)))
        % If the index corresponds to a variable:
        item.name = gep.variables(de_array(index) - length(gep.functions)); % Get the variable name.
        item.type = 2; % Type 2 indicates a variable.
        item.ports = 0; % Variables have no ports.
    else
        % If the index corresponds to a constant:
        item.name = gep.constants(de_array(index) - (length(gep.functions) + length(gep.variables))); % Get the constant value.
        item.type = 3; % Type 3 indicates a constant.
        item.ports = 0; % Constants have no ports.
    end
end
