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
## @deftypefn {Function File} ir_plotline( @var{beta}, @var{xlimits}, @var{ls} )
## Plots line y = beta(0) + beta(1) * x 
##
## @itemize @w
## @item Input
##   @itemize @w
##     @item @var{beta}:     (1 x 2)-matrix, line coefficients
##     @item @var{xlimits}:  vector of form [xmin xmax], limits for X axis
##     @item @var{ls}:       format string to control the line style.
##   @end itemize
## @end itemize
## @end deftypefn

function ir_plotline(beta, xlimits, ls)

if numel(beta) ~=2 
  error("Wrong dimension of beta");
endif

beta = beta(:);

x = xlimits(:);
X = [x.^0 x];
y = X * beta;

if nargin > 2
  plot(x, y, ls);
else
  plot(x, y);
end

end
