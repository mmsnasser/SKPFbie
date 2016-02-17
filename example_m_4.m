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
alphain =  0.0+0.2i;
%%
dv      = [0.5+0.5i;-0.5+0.5i;-0.5-0.5i; 0.5-0.5i];
qv      = [0.12    ; 0.11    ; 0.08    ; 0.15    ];
%%
t = [0:2*pi/n:2*pi-2*pi/n].';
delt  =  [0;dv];
q     =  [1;qv];
m     =   length(dv)
mn    =  (m+1)*n;
deltp =   delt./(abs(delt).^2-q.^2);
%%

%%
z       = [0.4+0.3i;-0.2-0.1i;+0.2-0.1i;-0.2+0.1i;0.8i];
alpha   =  0.0+1.90i;
%%
tic
vj;
toc
%%
tic
w   = skpf (z,alpha);
toc
%%
[z w]
%%
