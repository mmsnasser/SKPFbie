function  [mu , h ]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% The function 
%        [mu,h]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% return the unique solution mu of the integral equation 
%               (I-N)mu=-Mgam 
% and the function 
%                h=[(I-N)gam-Mmu]/2,
% where et is the parameterization of the boundary, etp=et', 
% A=exp(-i\thet)(et-alp) for bounded G and by A=exp(-i\thet) for unbounded
% G, gam is a given function, n is the number of nodes in  each boundary
% component, iprec is the FMM precision flag, restart is the maximum number 
% of GMRES method inner iterations, gmrestol is the tolerance of the GMRES   
% method, and maxit is the maximum number of GMRES method outer iterations
%
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

%%
a        = [real(et.') ; imag(et.')];
m        =  length(et)/n-1;
b1       = [etp./A].';
[Ub1]    = zfmm2dpart(iprec,(m+1)*n,a,b1,1);
Eone     = (Ub1.pot).';
%%
b(1,1) = 0;
for k=2:n
    b(k,1) = (-1)^(k+1)*(1/n)*cot(pi*(k-1)/n);
end
%%
[mu,~,~,~]      = gmres(@(x)fB(x),-fC(gam),restart,gmrestol,maxit);
if( nargout == 2 )
    h   = (fC(mu)-fB(gam))./2;
end
%%


%%
function  hx  = fB (x)
    bx2   = [x.*etp./A].';
    [Ubx2]= zfmm2dpart(iprec,(m+1)*n,a,bx2,1);
    Ex    = (Ubx2.pot).';
    hx    =  2.*x-(2/n).*imag(A.*Eone).*x+(2/n).*imag(A.*Ex);
end
%%
function  hx = fC (x)
    bx    = [x.*etp./A].';
    [Ubx] = zfmm2dpart(iprec,(m+1)*n,a,bx,1);
    Ex    = (Ubx.pot).';
    for k=1:m+1
        hLx(1+(k-1)*n:k*n,1) = ifft(fft(b).*fft(x(1+(k-1)*n:k*n,1)));
    end
    hx    = -(2/n).*real(A.*Ex)+(2/n).*real(A.*Eone).*x+hLx;    
end
%%
end