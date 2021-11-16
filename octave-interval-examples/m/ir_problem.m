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
## @deftypefn {Function File} {@var{irp} =} ir_problem(@var{X}, @var{y}, @var{epsilon}, @var{lb}, @var{ub}, @var{C}, @var{d}, @var{ctype})
##            {Function File} {@var{irp} =} ir_problem(@var{irpin}, @var{lb}, @var{ub}, @var{C}, @var{d}, @var{ctype})
##  Constructs interval regression problem for data and optional constraints
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{X}:  (n x m)-matrix of observations for X variables
##     @item @var{y}:  column n-vector of measured values of response variable
##     @item @var{epsilon}:  real number or n-vector, half-width(s) of uncertainty intervals for response variable @var{y}
##     @item @var{irpin}:  structure containing data for interval regression problem.
##     @item @var{lb}:  row m-vector containing the lower bound on each of the regression parameter. 
##                      If lb is not supplied, the default lower bound for the variables is -infinity.
##     @item @var{ub}:  row m-vector containing the upper bound on each of the regression parameter. 
##                      If lb is not supplied, the default upper bound for the parameters is assumed infinite.
##     @item @var{C}:  (k x m)-matrix of optional constraints coefficients for regression parameters beta: @var{C} * beta [ "=" |"<=" |">=" ] @var{d} 
##     @item @var{d}:  column k-vector of optional constraints right-hand side: @var{C} * beta [ "=" |"<=" |">=" ] @var{d}
##     @item @var{ctype}:  column k-vector of characters specifying optional constraints type:
##       @itemize @w
##         @item "U" An inequality constraint with an upper bound (@var{C(i,:)} * beta <= @var{d(i)}).
##         @item "S" An equality constraint (@var{C(i,:)} * beta = @var{d(i)}).
##         @item "L" An inequality with a lower bound (@var{C(i,:)} * beta >= @var{d(i)}).
##       @end itemize
##   @end itemize
## @item Output
##   @itemize @w
##     @item @var{irp}:  structure containing tha data and constrains for interval regression problem @var{y} = @var{X} * beta.
##   @end itemize
## @item Example
##   @itemize @w
##   @item X = [1 1; 1 1.5; 1 2];                    # X values
##   @item y = [1.5; 1.4; 2.5];                          # y values with bounded errors 
##   @item epsilon = 0.75;                           # y error bound
##   @item [irproblem] = ir_problem(X, y, epsilon);  # Define interval regression problem
##   @item figure
##   @item ir_scatter(irproblem);                    # Plot interval data
##   @item [beta, exitcode] = ir_outer(irproblem)    # Estimate regression parameters
##   @item beta =
##   @item 
##   @item -1.75000   2.75000
##   @item -0.50000   2.50000
##   @item 
##   @item exitcode =  1
##   @item 
##   @item [irproblem2] = ir_problem(irproblem, [-inf 0]); # Define interval regression problem with additional 
##   @item                                                 # lower bound constraints on regression parameters
##   @item 
##   @item [beta2, exitcode] = ir_outer(irproblem2)   # Estimate regression parameters
##   @item beta2 =
##   @item 
##   @item -1.75000   2.15000
##   @item  0.00000   2.50000
##   @item 
##   @item exitcode =  1
##   @end itemize
## @end itemize
##
## @end deftypefn

function [irp] = ir_problem(varargin)

if nargin < 1 
   error("Not enough arguments.");
   return;
endif

# Parse data for an interval regression problem 
if strcmp(class(varargin{1}), "struct")
   irp = varargin{1};
   varargin = varargin(2:end);
   nargin = nargin - 1;
else
  if nargin < 3
    error("Not enough arguments.");
    return;
  endif
  X = varargin{1};
  y = varargin{2};
  epsilon = varargin{3};

  if ~isIrpDimensionValid(X,y,epsilon)
    error("Wrong dimensions. Number of rows for X and y must be the same, dimension of epsilon must be equal to dimension of y or could be scalar.");
    return;
  endif
  
  if isempty(X) || isempty(y) || isempty(epsilon)
    error("X, y, epsilon must be non-empty");
    return;
  endif

  irp.X = X;
  irp.y = y;
  irp.epsilon = epsilon;

  varargin = varargin(4:end);
  nargin = nargin - 3;
endif

[n,m] = size(irp.X);

if nargin >= 1
  lb = varargin{1};
endif
if ~exist("lb","var") || isempty(lb)
  lb = -Inf * ones(1,m);
else
  if ~isempty(lb) && ~all(size(lb) == [1 m])
    error("Wrong dimensions of constraints: lb must be (1 x m)-matrix where m is the number of columns in X.");
    return;    
  endif
endif
irp.lb = lb;

if nargin >= 2
 ub = varargin{2};
endif
if ~exist("ub","var") || isempty(ub)
  ub = Inf * ones(1,m);;
else
  if ~isempty(ub) && ~all(size(ub) == [1 m])
    error("Wrong dimensions of constraints: ub must be (1 x m)-matrix where m is the number of columns in X.");
    return;
  endif  
endif
irp.ub = ub;

if nargin>=3 && nargin < 5
  error("Wrong number of parameters. Inequality constraints must be specified by C, d and ctype");
  return;
endif

if nargin >= 5
  C = varargin{3};
  d = varargin{4};
  ctype = varargin{5};
 
  if ~isempty(C) && size(C,2) ~= m 
    error("Wrong dimensions of constraints. Number of columns of C and X must be the same");
    return;
  endif

  if ~all(size(d) == size(C(:,1))) 
    error("Wrong dimensions of constraints. Dimensions of C(:,1) and d must be the same");
    return;
  endif

  if ~all(numel(ctype) == numel(d))
    error("Wrong dimensions of constraints. Dimensions of d and ctype must be the same");
    return;
  endif

  numUL = (upper(ctype) == 'L') + (upper(ctype) == 'U');
  if ~strcmp(class(ctype),"char") || sum(numUL(:)) ~= size(C,1)
    error("ctype must contain only 'U' or 'L' symbols");
    return;
  endif
else
  C = [];
  d = [];
  ctype = [];
endif

irp.C = C;
irp.d = d;
irp.ctype = ctype;

endfunction # ir_problem


function result = isIrpDimensionValid(X,y,epsilon)
  result = true;
 
  if size(y,1) ~= size(X,1) 
    result = false;
    return;
  endif
 
  if ~isscalar(epsilon) && ~all(size(epsilon) == size(y))
    result = false;
    return;    
  endif
 
endfunction
