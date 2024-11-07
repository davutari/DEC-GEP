% ReelLog(x) computes the natural logarithm of 
% the absolute value of the input x to avoid negative value x.
% It ensures real value results from log(x).
function [result]=ReelLog(x)
     result=log(abs(x));
end