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
## @deftypefn {Function File} ir_plotbeta( @var{irproblem}, @var{axisbox} )
## Plots 2D feasible parameters set for an interval regression problem
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}: structure with the description of interval 
##           regression problem 
##     @item @var{axisbox}: vector of form [xmin xmax ymin ymax] with axes limits 
##            (box is calculted automatically if it is omitted)   
##   @end itemize
## @item Example
##   @itemize @w
##     @item ## Simulate observations for simple dependency y = 1 + x 
##     @item n = 5;                                 # number of observations
##     @item epsilon = 0.5;                         # y error bound
##     @item X = [ ones(n,1) (1:5)' ];              # X values
##     @item y = X*[1 1]' + epsilon*rand(n,1)-0.5;  # y values with bounded errors 
##     @item [irproblem] = ir_problem(X, y, epsilon);
##     @item figure
##     @item ir_plotbeta(irproblem);
##   @end itemize
## @end itemize
##
## @end deftypefn

function ir_plotbeta(irproblem, axisbox)

X = irproblem.X;
y = irproblem.y;
epsilon = irproblem.epsilon;

if size(X,2) ~= 2
  error("Only 2D set of feasible parameters could be plotted")
endif
  
##  Unfold epsilon
if numel(epsilon)==1
   epsilon = epsilon*ones(numel(y),1);
endif
 
## Calculate vertecies of feasibile parameters set
vertices = ir_beta2poly(irproblem);
vertices = [vertices; vertices(1,:)];
 
actobs = [];

if nargin < 1 
  error("Too few arguments");
else
  if nargin == 1
    xymin = min(vertices);
    xymax = max(vertices);
    d = 0.1*(xymax-xymin);
    xymin = xymin - d;
    xymax = xymax + d;
    axisbox = [xymin(1) xymax(1) xymin(2) xymax(2)];
  endif
endif

XL=zeros(size(X,1),2);
YL=XL;
XH=XL;
YH=XL;

for i=1:size(X,1)
  if X(i,1)==0
    if X(i,2)==0
       disp("Zero x observation is not displayed");
       continue;
    endif
    XL(i,1:2) = axisbox(1:2);
    YL(i,:) = (y(i) - epsilon(i))/X(i,2);    
    XH(i,1:2) = axisbox(1:2);
    YH(i,:) = (y(i) + epsilon(i))/X(i,2);
  else
    XL(i,1) = (y(i) - epsilon(i) - X(i,2).*axisbox(3))./X(i,1);
    YL(i,1) = axisbox(3); 
    XH(i,1) = (y(i) + epsilon(i) - X(i,2).*axisbox(3))./X(i,1);
    YH(i,1) = axisbox(3); 
  endif
endfor

for i=1:size(X,1)
  if X(i,2)==0
    if X(i,1)==0
      continue;
    endif
    XL(i,:) = (y(i) - epsilon(i))/X(i,1); 
    YL(i,1:2) = axisbox(3:4);  
    XH(i,:) = (y(i) + epsilon(i))/X(i,1);
    YH(i,1:2) = axisbox(3:4);  
  else
    XL(i,2) = axisbox(1);
    YL(i,2) = (y(i) - epsilon(i) - X(i,1).*axisbox(1))./X(i,2);   
    XH(i,2) = axisbox(1);  
    YH(i,2) = (y(i) + epsilon(i) - X(i,1).*axisbox(1))./X(i,2);   
  endif
endfor

nactobs = setdiff(1:size(X,1),actobs);

hold on 

## Show non-active constraints
# plot(XL(nactobs,:)',YL(nactobs,:)',"k-");
# plot(XH(nactobs,:)',YH(nactobs,:)',"k-"); 

## Show active constraints
# plot(XL(actobs,:)',YL(actobs,:)',"r-","LineWidth",2.5);
# plot(XH(actobs,:)',YH(actobs,:)',"r-","LineWidth",2.5); 

## Show bounds of feasible parameters set 
patch(vertices(:,1),vertices(:,2),[0.72 0.72 0.9])
plot(vertices(:,1),vertices(:,2),"bo-","LineWidth",2,"MarkerSize",4,"MarkerFaceColor","b");

axis(axisbox);
box on

endfunction
