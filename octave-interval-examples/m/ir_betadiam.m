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
## @deftypefn {Function File} {[@var{diam}, @var{v1}, @var{v2}] =} ir_betadiam( @var{irproblem} )
## Calculates diameter of feasible prarameters set and end points of the diameter 
## (points from the set with maxiaml distance beteewn them)
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{irproblem}: structure with the description of interval
##           regression problem
##     @end itemize
## @item Output
##   @itemize @w
##     @item @var{diam}:      distance between points in feasible parameters set
##     @item @var{v1}, @var{v1}  n-vectors (where n is number of parameters),
##           points from feasible parameters set with maximal distance between 
##           them
##   @end itemize
## @end itemize
##
## @end deftypefn

function [diam, v1, v2] = ir_betadiam(irproblem)

vertices = ir_beta2poly(irproblem)';

D = sqdist(vertices,vertices);

[Dmax, idx] = max (D(:));
[r, c] = ind2sub (size (D), idx);

diam = sqrt(Dmax);
v1 = vertices(:,r)';
v2 = vertices(:,c)';

end


function D = sqdist(X1, X2)
% Pairwise square Euclidean distance between two sample sets
% Input:
%   X1, X2: dxn1 dxn2 sample matrices
% Output:
%   D: n1 x n2 square Euclidean distance matrix
% Written by Mo Chen (sth4nth@gmail.com).
D = bsxfun(@plus,dot(X2,X2,1),dot(X1,X1,1)')-2*(X1'*X2);

end
