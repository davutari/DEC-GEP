function [result] = div(x, y)
    epsilon = 0.001;
    % Checks if the denominator y is zero
    if y == 0
        % If y is zero, performs division by (y + epsilon) 
        % to avoid division by zero error.
        result = x ./ (y + epsilon);
    else
        % If y is not zero, perform normal division
        result = x ./ y;
    end
end

