function [U]=zfmm2dpart(iprec,nsource,source,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg)
%ZFMM2DPART Laplace particle target FMM in R^2 (Cauchy sums).
%
% Laplace FMM in R^2: evaluate all pairwise particle
% interactions (ignoring self-interactions) and interactions with targets.
%
% \phi(z_i) = \sum_{j\ne i}  dipstr_j *(1/(z_i-z_j))
%
% \gradient \phi(z_i) = \frac{\partial \phi(z_i)}{\partial z}
% \hessian  \phi(z_i) = \frac{\partial^2 \phi(z_i)}{\partial z^2}
%
% where dipstr are complex valued, z are complex numbers.
%
%
% [U]=ZFMM2DPART(IPREC,NSOURCE,SOURCE,DIPSTR);
%
% [U]=ZFMM2DPART(IPREC,NSOURCE,SOURCE,DIPSTR,...
%         IFPOT,IFGRAD,IFHESS);
%
% [U]=ZFMM2DPART(IPREC,NSOURCE,SOURCE,DIPSTR,...
%         IFPOT,IFGRAD,IFHESS,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG,IFHESSTARG);
%
%
% This subroutine evaluates the Cauchy potential, gradient and hessian due
% to a collection of dipoles. 
%
% Self-interactions are not-included.
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% nsource - number of sources
% source - real (2,nsource): source locations
% dipole - complex (nsource): dipole strengths
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% ifgrad - gradient computation flag, 1 => compute the gradient, otherwise no
% ifhess - hessian computation flag, 1 => compute the hessian, otherwise no
%
% ntarget - number of targets
% target - real (2,ntarget): target locations
%
% ifpottarg - target potential computation flag, 
%      1 => compute the target potential, otherwise no
% ifgradtarg - target gradient computation flag, 
%      1 => compute the target gradient, otherwise no
% ihesstarg - target hessian computation flag 
%      1 => compute the hessian, otherwise no
%
% Output parameters: 
%
% U.pot - complex (nsource) - potential at source locations
% U.grad - complex (2,nsource) - gradient  at source locations
% U.hess - complex (3,nsource) - hessian at source locations
% U.pottarg - complex (ntarget) - potential at target locations
% U.gradtarg - complex (2,ntarget) - gradient  at target locations
% U.hesstarg - complex (3,ntarget) - hessian at target locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%

if( nargin == 5 ) 
  ifpot = 1;
  ifgrad = 1;
  ifhess = 1;
  ntarget = 0;
  target = zeros(2,1);
  ifpottarg = 0;
  ifgradtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 7 ) 
  ntarget = 0;
  target = zeros(2,1);
  ifpottarg = 0;
  ifgradtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 9 ) 
  ifpottarg = 1;
  ifgradtarg = 1;
  ifhesstarg = 1;
end


ifpot = double(ifpot); ifgrad = double(ifgrad); ifhess = double(ifhess);
ifpottarg = double(ifpottarg); ifgradtarg = double(ifgradtarg);
ifhesstarg = double(ifhesstarg);

pot=0;
grad=0;
hess=0;
pottarg=0;
gradtarg=0;
hesstarg=0;

if( ifpot == 1 ), pot=zeros(1,nsource)+1i*zeros(1,nsource); end;
if( ifgrad == 1 ), grad=zeros(2,nsource)+1i*zeros(2,nsource); end;
if( ifhess == 1 ), hess=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifpottarg == 1 ), pottarg=zeros(1,ntarget)+1i*zeros(1,ntarget); end;
if( ifgradtarg == 1 ), gradtarg=zeros(2,ntarget)+1i*zeros(2,ntarget); end;
if( ifhesstarg == 1 ), hesstarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;

ier=0;

if( ntarget == 0 ) 
mex_id_ = 'zfmm2dpartself(io int[x], i int[x], i int[x], i double[xx], i dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pot, grad, hess] = fmm2d_r2012a(mex_id_, ier, iprec, nsource, source, dipstr, ifpot, pot, ifgrad, grad, ifhess, hess, 1, 1, 1, 2, nsource, 1, 1, 1);
else
mex_id_ = 'zfmm2dparttarg(io int[x], i int[x], i int[x], i double[xx], i dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pot, grad, hess, pottarg, gradtarg, hesstarg] = fmm2d_r2012a(mex_id_, ier, iprec, nsource, source, dipstr, ifpot, pot, ifgrad, grad, ifhess, hess, ntarget, target, ifpottarg, pottarg, ifgradtarg, gradtarg, ifhesstarg, hesstarg, 1, 1, 1, 2, nsource, 1, 1, 1, 1, 1, 1, 1);
end


if( ifpot == 1 ) U.pot=pot; end
if( ifgrad == 1 ) U.grad=grad; end
if( ifhess == 1 ) U.hess=hess; end
if( ifpottarg == 1 ) U.pottarg=pottarg; end
if( ifgradtarg == 1 ) U.gradtarg=gradtarg; end
if( ifhesstarg == 1 ) U.hesstarg=hesstarg; end
U.ier=ier;


