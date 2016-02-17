function  fz  = fcau (et,etp,f,z)
%%
% The function 
%        fz  = fcau (et,etp,f,z,n,finf)
% return the values of the analytic function f computed using the Cauchy
% integral formula at interior vector of points z, where et is the
% parameterization of the boundary, finf is the values of f at infinity 
% for unbounded G, n is the unber of nodes in each boundary component.
% The integral is discretized using the trapezoidal rule. The summations  
% are computed using the FMM.
% Copyright Mohamed Nasser, 2016
% Please cite this function as:
%  M.M.S. Nasser, Fast solution of boundary integral equations with the 
%  generalized Neumann kernel, Electronic Transactions on Numerical 
%  Analysis,  44 (2015) 189--229.
%
%
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.  This program is distributed in 
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details. You should have received a copy of the GNU General Public 
% License along with this program; 
% if not, see <http://www.gnu.org/licenses/>.
%%
vz    = [real(z) ; imag(z)];       % target
nz    = length(z);                 % ntarget
a     = [real(et.') ; imag(et.')]; % source
tn    = length(et);                % nsource=(m+1)n
iprec = 5;                         %- FMM precision flag
%%
bf    = [f.*etp].';
[Uf]  = zfmm2dpart(iprec,tn,a,bf,0,0,0,nz,vz,1,0,0);
b1    = [etp].';
[U1]  = zfmm2dpart(iprec,tn,a,b1,0,0,0,nz,vz,1,0,0);
%%
fz    = (Uf.pottarg)./(U1.pottarg);
%%
end