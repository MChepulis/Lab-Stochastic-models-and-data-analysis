## Copyright (C) 2021 Sergei Zhilin
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} ir_scatter( @var{irproblem}, @var{LineSpec} )
## Plots 2D interval scatter plot
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}: structure with the description of interval 
##           regression problem 
##     @item @var{LineSpec}: interval line style specification (default is "bo")
##   @end itemize
## @item Example
##   @itemize @w
##     @item ## Simulate observations for simple dependency y = 1 + x 
##     @item n = 5;                                  # number of observations
##     @item epsilon = 0.5;                          # y error bound
##     @item X = [ ones(n,1) (1:5)' ];               # X values
##     @item y = X*[1 1]' + epsilon*rand(n,1) - 0.5; # y values with bounded errors 
##     @item [irproblem] = ir_problem(X, y, epsilon);
##     @item figure
##     @item ir_scatter(irproblem,'ro');
##   @end itemize
## @end itemize
##
## @end deftypefn

function ir_scatter(irproblem, LineSpec)

if nargin < 2
    LineSpec="bo";
endif

X = irproblem.X;
y = irproblem.y;
epsilon = irproblem.epsilon;

if numel(epsilon) == 1
    epsilon = epsilon * ones(numel(y),1);
end

if size(X,2) == 2 && all(X(:,1)==1) 
    x = X(:,2);
elseif size(X,2) == 2 && all(X(:,2)==1) 
    x = X(:,1);
elseif size(X,1)==1 || size(X,2)==1
    x = X;
else
    error("Wrong X size");
end

errorbar(x,y,epsilon,LineSpec);

box on

end
