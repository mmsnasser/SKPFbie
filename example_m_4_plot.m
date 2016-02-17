clear, clear global
clc
%%
global n
global alphain
global dv
global qv
global mun
%%
n       = 2^9;
t       = [0:2*pi/n:2*pi-2*pi/n].';
alpha   =  0.0+0.4i;
alphain =  0.0+0.2i;
%%
dv      = [0.5+0.5i;-0.5+0.5i;-0.5-0.5i; 0.5-0.5i];
qv      = [0.12    ; 0.11    ; 0.08    ; 0.15    ];
%%
z  = [0.4+0.3i;-0.2-0.1i;+0.2-0.1i;-0.2+0.1i;0.8i];
%%
t = [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;dv];
q     =  [1;qv];
m     =   length(dv)
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
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
et                          =  1./conj(zet);
et(abs(et)>10)              =  NaN+i*NaN;
figure
hold on
for k=1:m+1
    c_cr    =  zet(1+(k-1)*n:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(c_cr,'k','LineWidth',2)
    c_cr    =  et(1+(k-1)*n:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(c_cr,':k','LineWidth',2)
end
plot(real(z),imag(z),'or','MarkerSize',2)
plot(real(alpha),imag(alpha),'db','MarkerFaceColor','b')
plot(real(alphain),imag(alphain),'sr','MarkerFaceColor','r')
axis equal
grid on
%%
aaa
%%
tic
mun  = vj;
omgz1  = skpf (z,alpha);
toc
%%

%%
tic
wf_L_7    = skprod(dv, qv, 6);
omgz_L_7  = wf_L_7(z, alpha);
toc
%%
format short g
abs(omgz_L_7-omgz1)
%%
