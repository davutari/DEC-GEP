function [line] = gep_line_create(size)
% Generates a line of specified size consisting of hyphens.
% Parameters:
%   - size: The desired length of the line.
% Returns:
%   - line: The generated line of hyphens.
  % Initialize an empty string to store the line
  line = [];
  % Loops through the specified size
  for i = 1:size
    % Appends a '-' character to the line string
    line = strcat(line, "-");
  end
end