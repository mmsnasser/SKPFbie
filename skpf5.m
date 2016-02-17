function  omgz  = skpf5 (mun ,delto , qo , n , z , alp , alpin , bound_test)
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
if abs(alp)>1
    alpst   = 1/conj(alp);
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
A        =  zet-alpin;
%%


%% 
for k=2:m+1
    if (abs(delt(k))>q(k))
        gam(:,k) = (1/(2*pi)).*log(abs((zet-delt(k))./(zet-deltp(k))));
    else
        gam(:,k) = (1/(2*pi)).*log(abs(zet-delt(k)));
    end        
end
for k=2:m+1
    if k==bound_test
        [~,hj(:,1)] = fbie(zet,zetp,A,gam(:,k),n,5,10,1e-14,10);
    end
end
%%

%%
walpst  = skpf4 (mun , delto , qo , n , z , alpst , alpin ,  bound_test);
%%


%%
j          =  bound_test;
fnj        = (gam(:,j)+hj+i.*mun(:,j))./A;
fnjz       =  fcau(zet,zetp,fnj,z.').';
%%
talpst     = -Arg((alpst-delt(bound_test))/q(bound_test),-2*pi);
munalpst   =  tripoly(mun(1+(j-1)*n:j*n,j),talpst);
%%


%%
if (abs(delt(j))>q(j))
    sHjaz    =  -1*(q(j)/(1-conj(delt(j))*alp)).*((z-delt(j))./(z-deltp(j)))...
           .*exp(-i*carg((alpst-delt(j))/(alpst-deltp(j)))+2*pi*i*munalpst...
           +2*pi*sum(hj(1:n))/n-2*pi.*(z-alpin).*fnjz);
else
    sHjaz    =  -1*(q(j)/(1-conj(delt(j))*alp)).*(z-delt(j))...
           .*exp(-i*carg(alpst-delt(j))+2*pi*i*munalpst...
           +2*pi*sum(hj(1:n))/n-2*pi.*(z-alpin).*fnjz);
end    
omgz    =  walpst./sHjaz;
%%
if (ztes==1 )
    if abs(alpo) == inf
        omgz = conj(omgz);
    elseif alpo ==0 
        omgz = conj(omgz);
    else
        omgz = -alpo.*zo.*conj(omgz);
    end
end
end