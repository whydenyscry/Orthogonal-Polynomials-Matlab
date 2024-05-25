function He = HermiteHePolynomialsSym(n, x)

if n == 0
    He = sym(1);
elseif n == 1
    He = [sym(1); x];
else
    He(1) = sym(1);
    He(2) = x;
    
    for i = 2:n
    He(i+1) = x.*He(i)-(i-1)*He(i-1);
    end
end

He = He';
end