function H = HermiteHPolynomialsSym(n, x)

if n == 0
    H = sym(1);
elseif n == 1
    H = [sym(1); 2 * x];
else
    H(1) = sym(1);
    H(2) = 2 * x;
    
    for i = 2:n
    H(i+1) = 2*x.*H(i)-2*(i-1)*H(i-1);
    end
end

H = H';
end