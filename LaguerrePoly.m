%-----------------------------------
%Laguerre Polynomial
%-----------------------------------

function y=LaguerrePoly(params,x)

n=params(1);
k=params(2);

m=[0:n];

a=factorial(n+k)*ones(1,length(m));
b=factorial(n-m);
c=factorial(k+m);
d=factorial(m);
e=(-1).^m;

y=zeros(size(x));
for s=1:n+1;
    y = y + a(s) ./ b(s) ./ c(s) ./ d(s) .* e(s) .* x.^m(s);
end
