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
## @deftypefn {Function File} ir_plotrect( @var{beta}, @var{ls} )
## Plots rectangle representing 2D outer interval estimate of feasible parameters 
## set for linear  regression y = X*beta with interval uncertainty in response 
## variable.
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{beta}:  (2 x 2)-matrix, lower (@var{beta(:,1)}) and upper 
##           bounds (@var{beta(:,2)}) for regression parameters estimates; 
##     @item @var{ls}:    format string to control the line style.
##   @end itemize
## @end itemize
## @end deftypefn

function ir_plotrect(beta,ls)
## This fuction is a shortcut for standard call:
## rectangle ("Position", [beta(:,1); beta(:,2)-beta(:,1)], 'EdgeColor', [1 0 0]);

if numel(beta) ~= 4
  error('beta must be 2 x 2 matrix')
endif

xidx = [1 1 3 3 1];
yidx = [2 4 4 2 2];
x = beta(xidx);
y = beta(yidx);

if nargin > 1
  plot(x,y,ls);
else
  plot(x,y);
endif

endfunction