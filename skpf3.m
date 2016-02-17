function  omgz  = skpf3 (mun , delto , qo , n , z , alp , alpin , bound_test)
%%
%
%
%%
ztes = 0;
if (abs(z)>1 )
    ztes = 1;
    zo   = z;
    alpo = alp;
    alp  = 1/conj(alp);
    z    = 1./conj(z); 
end
%%
t     = [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;delto];
q     =  [1;qo];
m     =   length(delto);
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
%%
zet(1:n,1)  = exp(1i.*t);
zetp(1:n,1) = 1i.*exp(1i.*t);
%%
for k=2:m+1
    zet(1+(k-1)*n:k*n,1)    = delt(k)+q(k).*exp(-1i.*t);
    zetp(1+(k-1)*n:k*n,1)   = -1i.*q(k).*exp(-1i.*t);
end
%%


%%
alpn     =  alpin;
A        =  zet-alpn;
%%

%%
if bound_test>0
    %%
    k = 1;
    phi(1+(k-1)*n:k*n,1)     =  0;
    for k=2:m+1
        phi(1+(k-1)*n:k*n,1) =  carg((zet(1+(k-1)*n:k*n)-deltp(k))./((zet(1+(k-1)*n:k*n)-alp).^2))...
                                +2*pi.*mun(1+(k-1)*n:k*n,k);
    end    
    %%
end
%%


%%
[psi , nu] =  fbie(zet,zetp,A,phi,n,5,10,1e-15,10);
fn         = (phi+nu+i.*psi)./A;
fnz        =  fcau(zet,zetp,fn,z.').';
%%

%%
ta      =  Arg((alp-delt(bound_test))/q(bound_test),0); %-2*pi
phialp  =  tripoly(phi(1+(bound_test-1)*n:(bound_test)*n,1),ta);
psialp  =  tripoly(psi(1+(bound_test-1)*n:(bound_test)*n,1),ta);
nualp   =  sum(     nu(1+(bound_test-1)*n:(bound_test)*n,1))/n;
Aalp    =  alp-alpin;
fnalp   = (phialp+nualp+i.*psialp)/Aalp;
const_C =  exp(-1i*(alp-alpin).*fnalp./2);
%%

%%
omgz    =  const_C.*(z-alp).*exp(1i*(z-alpin).*fnz./2);
%%

%%
if (ztes==1 )
    omgz = -alpo.*zo.*conj(omgz);
end
end