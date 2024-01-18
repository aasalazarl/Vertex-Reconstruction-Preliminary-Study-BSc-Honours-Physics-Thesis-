function [c,norm_r] = linefit3D(x, y, z)
% linefit    Least-squares fit of data to z = c(1)*x + c(2)*y + c*(3)
%
% Synopsis:   c     = linefit(x,y)
%            [c,norm_r] = linefit(x,y)
%
% Input:   x,y = vectors of independent and dependent variables
%
% Output:  c  = vector of slope, c(1), and intercept, c(2) of 
%               least sq. line fit
%          norm_r = (optional) norm of residual


if length(z) == length(x) && length(z) == length(y)
    x = x(:);  y = y(:); z = z(:);    %  Make sure that x and y 
                                      % are column vectors
    A = [x y ones(size(x))];  %  m-by-n matrix of overdetermined system
    c = (A'*A)\(A'*z);      %  Solve normal equations
    
    %if nargout>1
     %   r = y - A*c;
      %  norm_r = norm(r);
    %end

else
    error('x and y are not compatible'); 
end