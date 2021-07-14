%-----------------------------------
% Hermite Polynomial
%-----------------------------------

function y=hermitePoly(n,x)

m=(0:floor(n/2));

a=factorial(n-2*m);
b=factorial(m);

y=zeros(size(x));
for s=1:length(m)
    y = y   +   factorial(n) ./ a(s) ./ b(s)  .*   (-1).^m(s)  *  (2*x).^(n-2*m(s));
end
