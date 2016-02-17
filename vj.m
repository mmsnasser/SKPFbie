function  mun   = vj 
%%
%
%
%%
global n
global alphain
global dv
global qv
global mun
%%
t     =  [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;dv];
q     =  [0;qv];
m     =   length(dv);
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
%%
zet(1:n,1)  = exp(1i.*t);
zetp(1:n,1) = 1i.*exp(1i.*t);    
for k=2:m+1
    zet(1+(k-1)*n:k*n,1)    = delt(k)+q(k).*exp(-1i.*t);
    zetp(1+(k-1)*n:k*n,1)   = -1i.*q(k).*exp(-1i.*t);
end
%%
A        =  zet-alphain;
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
    mun(:,k) = fbie(zet,zetp,A,gam(:,k),n,5,10,1e-14,10);
end
%%
end