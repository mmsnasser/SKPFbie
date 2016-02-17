clear, clear global
clc
%%
n       = 2^10;
t       = [0:2*pi/n:2*pi-2*pi/n].';
alpha   =  0.1+0.6i;
alphain =  0.1+0.1i;
%%
dv      = [-0.76-0.28i;-0.76-0.04i;-0.76+0.20i;-0.52-0.52i;-0.52-0.28i;...
           -0.52-0.04i;-0.52+0.20i;-0.52+0.44i;-0.52+0.68i;-0.28-0.76i;...
           -0.28-0.52i;-0.28-0.28i;-0.28-0.04i;-0.28+0.20i;-0.28+0.44i;...
           -0.28+0.68i;-0.04-0.76i;-0.04-0.52i;-0.04-0.28i;-0.04-0.04i;...
           -0.04+0.20i;-0.04+0.44i;-0.04+0.68i; 0.20-0.76i; 0.20-0.52i;...
            0.20-0.28i; 0.20-0.04i; 0.20+0.20i; 0.20+0.44i; 0.20+0.68i;...
            0.44-0.52i; 0.44-0.28i; 0.44-0.04i; 0.44+0.20i; 0.44+0.44i;...
            0.44+0.68i; 0.68-0.52i; 0.68-0.28i; 0.68-0.04i; 0.68+0.20i;...
            0.68+0.44i];

qv      = [ 0.082502  ; 0.075771  ; 0.093485  ; 0.091531  ; 0.09798   ;...
            0.055054  ; 0.09025   ; 0.07843   ; 0.064497  ; 0.035358  ;...
            0.097875  ; 0.05937   ; 0.082884  ; 0.094489  ; 0.039119  ;...
            0.041126  ; 0.07781   ; 0.072215  ; 0.10747   ; 0.09951   ;...
            0.093888  ; 0.034466  ; 0.036341  ; 0.037702  ; 0.099457  ;...
            0.11204   ; 0.089483  ; 0.041491  ; 0.092877  ; 0.039601  ;...
            0.040222  ; 0.085742  ; 0.058607  ; 0.086882  ; 0.095174  ;...
            0.080737  ; 0.094383  ; 0.05043   ; 0.093941  ; 0.11444   ;...
            0.10542   ];
%%



%%
t = [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;dv];
q     =  [1;qv];
m     =   length(dv)
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
%%
[xzo,yzo] = meshgrid([-1.0:0.2:1],[-1.0:0.2:1]);
zo = xzo+1i*yzo;
zo(abs(zo)>=1)=NaN+1i*NaN;
for k=1:m
    zo(abs(zo-dv(k))<=qv(k))=NaN+1i*NaN;
end
z1 =  zo(abs(zo)>=0).';
z2 =  1./conj(z1);
z  = [z1,z2];
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
zeto                          =  1./conj(zet);
zeto(abs(zeto)>50)            =  NaN+i*NaN;
%%

%%
figure
hold on
for k=1:m+1
    c_cr    =  zet(1+(k-1)*n:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(c_cr,'k','LineWidth',2)
end
plot(real(z),imag(z),'or','MarkerSize',2)
plot(real(alpha),imag(alpha),'db','MarkerFaceColor','b')
plot(real(alphain),imag(alphain),'sr','MarkerFaceColor','r')
axis equal
grid on
%%
figure
hold on
for k=1:m+1
    c_cr    =  zet(1+(k-1)*n:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(c_cr,'k','LineWidth',2)
    c_cr    =  zeto(1+(k-1)*n:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(c_cr,':k','LineWidth',2)
end
plot(real(z),imag(z),'or','MarkerSize',2)
plot(real(alpha),imag(alpha),'db','MarkerFaceColor','b')
plot(real(alphain),imag(alphain),'sr','MarkerFaceColor','r')
axis equal
grid on
%%