function [y,yp,ypp] = tripoly (x,s)
% Suppose that x is an 1*n vector. The vector x is the values of f(t)
% at the n equidistant points
% t_j=(j-1)*2*pi/n,  j=1,2,...,n.
% This function find the trigonometric interpolating polynomial that
% interpolates the function f(t) at the n points t_1,t_2,...,t_n.
% The function compute the values of the function f(t) at the point s
% (s can be an m*1 vector of points), its first derivative, and its
% second derivative, i.e.
% y = f(s),  yp=f'(s),  ypp=f''(s)
%%
n  = max(size(x));
A = fft(x)./n;
%%
a(1)=A(1);
b(1)=0;
for k=2:n/2+1;
    a(k)=2*real(A(k));
    b(k)=-2.*imag(A(k));
end
a(n/2+1)=0.5*a(n/2+1);
b(n/2+1)=0;
%%
y = a(1)+s-s;
for k=1:n/2
    y  = y+(a(k+1).*cos(k.*s)+b(k+1).*sin(k.*s));
end
%%
if nargout>1
yp = s-s;
for k=1:n/2
    yp  = yp+(-k*a(k+1).*sin(k.*s)+k*b(k+1).*cos(k.*s));
end
end
%%
if nargout==3
ypp = s-s;
for k=1:n/2
    ypp  = ypp+(-(k^2)*a(k+1).*cos(k.*s)-(k^2)*b(k+1).*sin(k.*s));
end
end
%%
end