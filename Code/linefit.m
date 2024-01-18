function [c,norm_r] = linefit(x,y)
% linefit    Least-squares fit of data to y = c(1)*x + c(2)
%
% Synopsis:   c     = linefit(x,y)
%            [c,norm_r] = linefit(x,y)
%
% Input:   x,y = vectors of independent and dependent variables
%
% Output:  c  = vector of slope, c(1), and intercept, c(2)
%               of least sq. line fit
%          norm_r = (optional) norm of residual

if length(y)~= length(x),  error('x and y are not compatible');  end

x = x(:);  y = y(:);    %  Make sure that x and y are column vectors
A = [x ones(size(x))];  %  m-by-n matrix of overdetermined system
c = (A'*A)\(A'*y);      %  Solve normal equations
if nargout>1
  r = y - A*c;
  norm_r = norm(r);
end
