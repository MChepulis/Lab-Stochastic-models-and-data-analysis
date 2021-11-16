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
## @deftypefn {Function File} {[@var{beta}, @var{exitcode}, @var{active}] =} ir_outer( @var{irproblem} )
## Builds outer box estimate of feasible parameters set for linear 
## regression y = X*beta with interval uncertainty in response variable
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}: structure with the description of interval regression problem containing 
##           data (irproblem.X, irproblem.y, irproblem.epsilon) and optional constrains for beta 
##           (irproblem.C, irproblem.d, irproblem.ctype, irproblem.ul, irproblem.ub).
##   @end itemize
## @item Output
##   @itemize @w
##     @item @var{beta}:  (m x 2)-matrix, lower (@var{beta(:,1)}) and upper bounds (@var{beta(:,2)})
##           for regression parameters estimates; m is number of columns in matrix @var{irproblem.X}
##     @item @var{exitcode} integer exit code
##        @itemize @w
##           @item  1:  Ok, solution is found
##           @item -1:  Unbounded solution set (colinearity in the data)
##           @item -2:  No feasible solution (feasible parameters set is empty)
##        @end itemize
##     @item @var{active}: struct with fields @var{active.lower} and @var{active.upper} 
##           containing vectors of active constraints (observations) which form 
##           the boundary of feasible parameters set
##   @end itemize
## @item Example
##   @itemize @w
##   @item ## Simulate observations for simple dependency y = 1 + x 
##   @item n = 5;                                 # number of observations
##   @item epsilon = 0.5;                         # y error bound
##   @item X = [ ones(n,1) (1:5)' ];              # X values
##   @item y = X*[1 1]' + epsilon*rand(n,1)-0.5;  # y values with bounded errors 
##   @item [irproblem] = ir_problem(X, y, epsilon);
##   @item figure
##   @item ir_scatter(irproblem);
##   @item [beta, exitcode, active] = ir_outer(irproblem) # estimate regression parameters
## @end itemize
##
## @end deftypefn

function [beta, exitcode, active] = ir_outer(irproblem)

X = irproblem.X;
y = irproblem.y;
epsilon = irproblem.epsilon;

C = irproblem.C;
d = irproblem.d;
Cdctype = irproblem.ctype;

if ~isfield(irproblem,'lb') || isempty(irproblem.lb)
  lb = -Inf*ones(size(X,2),1);
else
  lb = irproblem.lb;
endif

if ~isfield(irproblem,'ub') || isempty(irproblem.ub)
  ub = Inf*ones(size(X,2),1);
else
  ub = irproblem.ub;
endif

## Build matrix and right-hand side vector of linear programming problem
A = [X; -X; C];
b = [y+epsilon; -y+epsilon; d];
[n m] = size(X);

ctype(1:2*n) = 'U'; # constraints type is inequality
if ~isempty(Cdctype)
  ctype = [ctype; Cdctype];
endif
vartype(1:m) = 'C'; # variable type is continuous
sense = 1; % find minimum
  
## Siginificant digit for active bounds detection
SIGNIFICANT=0.0000001;
  
## Parameters estiamtes
beta = [];

## Lagrange multipliers (to detect active constraints)
L = [];

## Solve 2*m linear programming problems in order to find 
## lower and upper bounds for feasible parameters set on each axis
for i = 1:m
  
  ## Set objective function coefficients
  c = zeros(1,m);
  c(i) = 1;

  ## Solve LPPs and save indicies for non-zero Lagrange multipliers

  exitcode = 1;
  [btlow, flow, errcode, lambda] = glpk (c, A, b, lb, ub, ctype, vartype, sense);
  if errcode > 0 
      exitcode = -errcode;
      return
  endif
  L = unique([L; find(abs(lambda.lambda) > SIGNIFICANT)]);

  [bthigh, fhigh, errcode, lambda] = glpk (-c, A, b, lb, ub, ctype, vartype, sense);
  if errcode > 0 
      exitcode = -errcode;
      return
  endif
  L = unique([L; find(abs(lambda.lambda) > SIGNIFICANT)]);
  
  ## Save lower and upper bounds for i-th parameter
  beta = [beta; [flow -fhigh]];
endfor

## Detect active boundary observations
active.lower = L(find(L>n))-n;
active.upper = L(find(L<=n));

## TODO: add active optional constraints 

endfunction
