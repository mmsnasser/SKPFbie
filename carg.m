function  y   =   carg (z)
%This function compute the (continuous function) arg(f)
%
%
%%
n       =  length(z);
y       =  unwrap(angle(z));
%%
end