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
## @deftypefn {Function File} {[@var{yp}, @var{betap}, @var{exitcode}, @var{active}]} ir_predict(@var{irproblem}, @var{Xp})
##  Computes interval prediction for response @var{y} at the points @var{Xp} 
##  using the interval regression y = X*beta constructed for problem @var{irproblem}
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{Xp}:        (k x m)-matrix, contains points (row vectors) to predict at
##     @item @var{irproblem}: structure describing interval regression problem 
##   @end itemize
## @item Output
##   @itemize @w
##     @item @var{yp}:        (k x 2)-matrix, interval regression predictions 
##            at points Xp(i,:), yp(i,1) and yp(i,2) are, accordingly, lower 
##            and upper bounds of predicted interval   
##     @item @var{betap}:     (k x m x 2)-matrix, regression parameters values 
##            providing interval prediction yp(i), beta(i,:,1) and beta(i,:,2)
##            correspond to lower and upper bounds of yp accordingly
##     @item @var{exitcode}:  integer exit code
##        @itemize @w
##          @item  1  Ok, solution is found
##          @item -1  Unbounded solution set (colinearity in the data)
##          @item -2  No feasible solution (feasible parameters set is empty)
##        @end itemize
##     @item @var{active}     k-element array of structs with fields 'lower' 
##             and 'upper', active(i).lower and active(i).upper contain indicies 
##             of active  constraints (observations) which limit the lower and 
##             upper bounds of a predicted interval yp(i,:)
##   @end itemize
## @item Example
##   @itemize @w
##   @item X = [1 1; 1 1.5; 1 2];                    # X values
##   @item y = [1.5; 1.4; 2.5];                      # y values with bounded errors 
##   @item epsilon = 0.75;                           # y error bound
##   @item [irproblem] = ir_problem(X, y, epsilon);  # Define interval regression problem
##   @item figure
##   @item ir_scatter(X,y,epsilon);                  # Show interval dataset
##   @item Xp = [1 0.5; 1 1.75];                     # Set points where to predict
##   @item [yp, betap, exitcode] = ir_predict(irpproblem, Xp)
##   @item ypmid = mean(yp,2);                       # Get intervals' middle points
##   @item yprad = (yp(:,2)-yp(:,1))/2;              # Get intervals' radii
##   @item hold on
##   @item ir_scatter(x,ypmid,yprad,'b.');           # Show predicted intervals
##   @end itemize
## @end itemize
##
## @end deftypefn

function [yp, betap, exitcode, active] = ir_predict(irproblem, Xp)

SIGNIFICANT = 0.00000001;

X = irproblem.X;
y = irproblem.y;
epsilon = irproblem.epsilon;

lb = irproblem.lb;
ub = irproblem.ub;

C = irproblem.C;
d = irproblem.d;
Cdctype = irproblem.ctype;

if size(Xp,2) ~= size(X,2)
    error("Xp must have the same number of columns as X");
endif

A = [X; -X; C];
b = [y+epsilon; -y+epsilon; d];

k = size(Xp,1);
[n, m] = size(X);

ctype(1:2*n) = 'U'; # constraints type is inequality
if ~isempty(Cdctype)
  ctype = [ctype; Cdctype];
endif
vartype(1:m) = 'C'; # variable type is continuous
sense = 1;          # find minimum

## Allocate matricies and structures
yp = zeros(size(Xp,1),2);
betaplow = zeros(size(Xp));
betaphigh = betaplow;
active = struct('lower',{},'upper',{});

for i = 1:k
    [betalow, flow, exitcode, lambda] = glpk (Xp(i,:), A, b, lb, ub, ctype, vartype, sense);
    actlow=[];
    for j = 1:size(lambda.lambda,1)
        if abs(lambda.lambda(j)) > SIGNIFICANT
            actlow = [actlow j]; 
        endif
    endfor

    [betahigh, fhigh, exitcode, lambda] = glpk (-Xp(i,:), A, b, lb, ub, ctype, vartype, sense);
    actupp=[];
    for j = 1:size(lambda.lambda,1)
        if abs(lambda.lambda(j)) > SIGNIFICANT
            actupp = [actupp j]; 
        endif
    endfor

    yp(i,:) = [flow, -fhigh];
    betaplow(i,:)  = min(betalow,betahigh);
    betaphigh(i,:) = max(betalow,betahigh);
    active(i).lower = actlow;
    active(i).upper = actupp;
endfor

betap = cat(3, betaplow, betaphigh);

endfunction 