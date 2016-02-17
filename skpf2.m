function  omgb  = skpf2 (mun , delto , qo , n , alp , alpin , near_test)
%%
%
%
%%
t = [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;delto];
q     =  [0;qo];
m     =   length(delto);
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
%%
thet_fun  =  @(z,j)(delt(j)+q(j)^2.*z./(1-conj(delt(j)).*z));
%%

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
alpst    = 1/conj(alp);
A        =  zet-alpin;
%%

%% 
if (alp==0 )
    gam(:,1)  = (1/(2*pi)).*log(abs(zet));
elseif (abs(alp)==inf)
    gam(:,1)  = -(1/(2*pi)).*log(abs(zet));
elseif near_test<=1
    gam(:,1)  = (1/(2*pi)).*log(abs((zet-alp)./(zet-alpst)));
else
    gam(:,1)  = (1/(2*pi)).*log(abs((zet-alp).*(zet-thet_fun(alp,near_test))./((zet-alpst).*(zet-thet_fun(alpst,near_test)))));
end
mun(:,1)      =  fbie(zet,zetp,A,gam,n,5,10,1e-14,10);
%%

%%            
if near_test<=1
    phi   = -2*pi*mun(:,1);
    k = 1;
    if (alp==0 | abs(alp)==inf)
        phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+0;
    else
        phi(1+(k-1)*n:k*n,1)     =  phi(1+(k-1)*n:k*n,1)+...
                                carg(alp.*zet(1+(k-1)*n:k*n)./((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst)));
    end
    for k=2:m+1
        if (abs(delt(k))>q(k))
            if (alp==0 | abs(alp)==inf)
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg((zet(1+(k-1)*n:k*n)-deltp(k))./(zet(1+(k-1)*n:k*n)));        
            else
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg(alp.*(zet(1+(k-1)*n:k*n)-deltp(k))./((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst)));
            end
        else
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg(alp./((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst)));            
        end            
    end   
else
    phi   = -2*pi*mun(:,1);
    k = 1;
    if (alp==0 | abs(alp)==inf)
        phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+...
                                carg((zet(1+(k-1)*n:k*n)-thet_fun(zet(1+(k-1)*n:k*n),near_test)).^2./...
                                ((zet(1+(k-1)*n:k*n)-thet_fun(alp,near_test)).*(zet(1+(k-1)*n:k*n)-thet_fun(alpst,near_test))));
    else
        phi(1+(k-1)*n:k*n,1)     =  phi(1+(k-1)*n:k*n,1)+...
                                carg(alp.*zet(1+(k-1)*n:k*n).*(zet(1+(k-1)*n:k*n)-thet_fun(zet(1+(k-1)*n:k*n),near_test)).^2./...
                                ((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst).*...
                                (zet(1+(k-1)*n:k*n)-thet_fun(alp,near_test)).*(zet(1+(k-1)*n:k*n)-thet_fun(alpst,near_test))));
    end
    for k=2:m+1
        if (abs(delt(k))>q(k))
            if (alp==0 | abs(alp)==inf)
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg((zet(1+(k-1)*n:k*n)-deltp(k)).*(zet(1+(k-1)*n:k*n)-thet_fun(zet(1+(k-1)*n:k*n),near_test)).^2./...
                                (zet(1+(k-1)*n:k*n).*(zet(1+(k-1)*n:k*n)-thet_fun(alp,near_test)).*(zet(1+(k-1)*n:k*n)-thet_fun(alpst,near_test))));        
            else
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg(alp.*(zet(1+(k-1)*n:k*n)-deltp(k)).*(zet(1+(k-1)*n:k*n)-thet_fun(zet(1+(k-1)*n:k*n),near_test)).^2./...
                                ((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst).*...
                                (zet(1+(k-1)*n:k*n)-thet_fun(alp,near_test)).*(zet(1+(k-1)*n:k*n)-thet_fun(alpst,near_test))));
            end
        else
                phi(1+(k-1)*n:k*n,1) =  phi(1+(k-1)*n:k*n,1)+2*pi.*mun(1+(k-1)*n:k*n,k)+...
                                carg(alp.*(zet(1+(k-1)*n:k*n)-thet_fun(zet(1+(k-1)*n:k*n),near_test)).^2./...
                                ((zet(1+(k-1)*n:k*n)-alp).*(zet(1+(k-1)*n:k*n)-alpst).*...
                                (zet(1+(k-1)*n:k*n)-thet_fun(alp,near_test)).*(zet(1+(k-1)*n:k*n)-thet_fun(alpst,near_test))));
        end            
    end
end
%%

%%
[psi , nu] =  fbie(zet,zetp,A,phi,n,5,10,1e-15,10);
fn         = (phi+nu+i.*psi)./A;
fnalp      =  fcau(zet,zetp,fn,alp);
const_C    =  exp(-1i*(alp-alpin).*fnalp./2);
if near_test<=1
    omgb    =  const_C.*(zet-alp).*exp(1i*(zet-alpin).*fn./2);
else
    omgb    =  const_C.*(zet-alp).*((zet-thet_fun(alp,near_test))./(zet-thet_fun(zet,near_test))).*exp(1i*(zet-alpin).*fn./2);
end
%%
end