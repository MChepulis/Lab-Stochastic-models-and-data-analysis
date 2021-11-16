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
## @deftypefn {Function File} ir_plotmodelset( @var{irproblem}, @var{xlimits} )
## Plots 2D set of feasible linear models consistent with dataset of IR problem
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}: structure with the description of interval 
##           regression problem 
##     @item @var{xlimits}: vector of form [xmin xmax], limits for X axis
##   @end itemize
## @item Example
##   @itemize @w
##   @item ## Define an interval regression problem
##   @item X = [1 1; 1 1.5; 1 2];                   # X values
##   @item y = [1.5; 1.4; 2.5];                     # y values 
##   @item epsilon = 0.75;                          # y error bound
##   @item [irproblem] = ir_problem(X, y, epsilon); # Define IR problem
##   @item figure
##   @item ir_plotmodelset(irproblem);              # Plot set of models
##   @item
##   @item ## Define another interval regression problem with additional 
##   @item ## lower bound constraints on regression parameters
##   @item [irproblem2] = ir_problem(irproblem, [-inf 0]); 
##   @item 
##   @item figure
##   @item ir_plotmodelset(irproblem2);             # Plot another set of models
##   @item 
##   @end itemize
## @end itemize
##
## @end deftypefn

function ir_plotmodelset(irproblem, xlimits)
   X = irproblem.X;
   y = irproblem.y;
   epsilon = irproblem.epsilon;

   if size(X,1) < 2
      error("Not enough data");
   endif

   if size(X,2) == 2 && all(X(:,1)==1) 
       onescol = 1;
       xcol = 2;
   elseif size(X,2) == 2 && all(X(:,2)==1) 
       xcol = 1;
       onescol = 2;
   elseif size(X,2)==1
       xcol = 1;
   else
       error("Wrong X size");
   endif

   if ~exist("xlimits","var")
      Xbefore = 2*X(1,:) - X(2,:);
      Xafter  = 2*X(end,:) - X(end-1,:);
   else
      Xbefore = X(1,:);
      Xbefore(xcol) = xlimits(1);
      Xafter = X(1,:);
      Xafter(xcol) = xlimits(2);
   endif
   Xp = [Xbefore; X; Xafter];

   x = Xp(:,xcol);
   [yp, betap, exitcode] = ir_predict(irproblem, Xp);
   px = [x; flipud(x)];
   py = [yp(:,1); flipud(yp(:,2))];
   pcolor = [250 197 250]/255;
   # patch(px,py,pcolor);
   hold on
   plot(x,yp(:,1),"m-","LineWidth",1);
   plot(x,yp(:,2),"m-","LineWidth",1);
endfunction
   