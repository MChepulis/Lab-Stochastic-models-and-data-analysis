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
## @deftypefn {Function File} {[@var{vertices}] =} ir_beta2poly( @var{irproblem} )
## Calculates vertices of feasible parameters set polytope 
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}:  structure with the description of interval
##           regression problem
##     @end itemize
## @item Output
##   @itemize @w
##     @item @var{vertices}:   (m x n)-matrix, m row vectors representing 
##           n-dimensional vertices of feasible parameters set polytope
##   @end itemize
## @end itemize
##
## @end deftypefn

function [vertices] = ir_beta2poly(irproblem)

## Get IR problem data
X = irproblem.X;
y = irproblem.y;
epsilon = irproblem.epsilon;
C = irproblem.C;
d = irproblem.d;
ctype = irproblem.ctype;
lb = irproblem.lb;
ub = irproblem.ub;

##  Unfold epsilon
if numel(epsilon)==1
  epsilon = epsilon*ones(numel(y),1);
endif
 
## Prepare inequalities defining feasible parameters set polytope
A = [X; -X];
b = [y+epsilon; -y+epsilon];
 
Lidx = (ctype == 'L');
A = [ A; -C(Lidx,:) ]; 
b = [ b; -d(Lidx) ];

Sidx = (ctype == 'S');
Aeq = C(Sidx,:); 
beq = d(Sidx);

[A,b,Aeq,beq] = addBounds(A,b,Aeq,beq,lb,ub);
 
## Calculate vertices of feasibile parameters set
warning('off', 'Octave:colon-nonscalar-argument');
[V,nr] = lcon2vert(A,b,Aeq,beq); 

dim = size(V,2);
if dim == 1
  k = nr;
elseif dim == 2  
  k = convhull(V(:,1),V(:,2));
  k = k(1:end-1);
else
  h = convhulln(V);
  k = unique(h(:));   
endif

vertices = V(k,:);
 
endfunction
