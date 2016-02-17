function  omgz   = skpf  (z , alp )
%
%
%
%
% skpf function for computing the S-K prime function using a boundary
% integral method based on the boundary integral equation with the
% generalized Neumann kernel.
%
% -----
% Input:
%   dv = a vector of circle centers.
%   qq = a vector of circle radii.
%   n  =   the number of node in the discretization of each boundary circle.
%   alphain = a select point in the bounded domain D.
%   dv, dq, n, alphain: must be global
%
%   Then, we call:
%   vj;
%   omgz   = skpf (z,alpha);
%
% Output:
%     w = wf(z, alpha),
%   where z is an array of points in the bounded domain D or in the 
%   unbouded domain D' (not on the boundary) at which to evaluate the 
%   function, and alpha is a scalar parameter value (can be any point
%   except the two cases: alpha=inf and lies on the boundary or alpha=0 
%   and lies on the boundary.
%
%
% Copyright Mohamed Nasser, 2016
% Please cite this function as:
% M. M. S. Nasser, SKPFbie: a fast boundary integral equation 
% implementation of the Schottky–Klein prime function in MATLAB, 
% Version 1.1, 2016. https://github.com/mmsnasser/SKPFbie.
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
global n
global alphain
global dv
global qv
global mun
%%
m     =   length(dv);
delt  =  [0;dv];
q     =  [1;qv];
bound_test = 0;
if (abs(alp)<=1)
    for k=1:m+1
        if abs(abs(alp-delt(k))-q(k))<1e-15
            bound_test = k;
        end
    end
end
if (abs(alp)>1)
    alpst   = 1/conj(alp);
    for k=1:m+1
        if abs(abs(alpst-delt(k))-q(k))<1e-15
            bound_test = k;
        end
    end
end
%%
near_test  = 0;
if (abs(alp)<1)
    for k=1:m+1
        if abs(abs(alp-delt(k))-q(k))>1e-15 & abs(abs(alp-delt(k))-q(k))<=1e-2
            near_test = k;
        end
    end
end
if (abs(alp)>1)
    alpst   = 1/conj(alp);
    for k=1:m+1
        if abs(abs(alpst-delt(k))-q(k))>1e-15 & abs(abs(alpst-delt(k))-q(k))<=1e-2
            near_test = k;
        end
    end
end
%%

%%
ki   =  0; ko   =  0; kb  =  0;
for k=1:length(z)
    if abs(z(k))<1
        ki = ki+1;
        zi(ki,1)=z(k);
    end
    if abs(z(k))>1
        ko = ko+1;
        zo(ko,1)=z(k);
    end
    if abs(z(k))==1
        'Error: z must be interior point.'
        'Next version will consider boundary points.'
    end
end
%%
alpha_on_bundary  =  0;
if abs(abs(alp)-1)<1e-15
    if ki>0
        omgzi   = skpf3 (mun,dv,qv,n,zi,alp,alphain,bound_test);
    end
    if ko>0
        omgzo   = skpf3 (mun,dv,qv,n,zo,alp,alphain,bound_test);
    end
    alpha_on_bundary = 1;
end

for k=1:length(qv)
    alpst   = 1/conj(alp);
    if abs(abs(alp-dv(k))-qv(k))<1e-15
        if ki>0
            omgzi   = skpf4 (mun,dv,qv,n,zi,alp,alphain,bound_test);
        end
        if ko>0
            omgzo   = skpf5 (mun,dv,qv,n,zo,alp,alphain,bound_test);
        end
        alpha_on_bundary = 1;
    end
    if abs(abs(alpst-dv(k))-qv(k))<1e-15
        if ki>0
            omgzi   = skpf5 (mun,dv,qv,n,zi,alp,alphain,bound_test);
        end
        if ko>0
            omgzo   = skpf4 (mun,dv,qv,n,zo,alp,alphain,bound_test);
        end
        alpha_on_bundary = 1;
    end
end

if alpha_on_bundary ==0
    if ki>0
        omgzi   = skpf1  (mun,dv , qv , n , zi , alp, alphain, near_test);
    end
    if ko>0
        omgzo   = skpf1  (mun,dv , qv , n , zo , alp, alphain, near_test);
    end
    if ki==0 & ko==0
        'what is z?'
    end
end
%%
ki   =  0; ko   =  0; kb  =  0;
for k=1:length(z)
    if abs(z(k))<1
        ki = ki+1;
        omgz(k,1) = omgzi(ki);
    end
    if abs(z(k))>1
        ko = ko+1;
        omgz(k,1) = omgzo(ko);
    end
end
%%
end