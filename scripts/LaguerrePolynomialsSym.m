function L = LaguerrePolynomialsSym(alpha, n, x)

if n == 0
    L = 1;
elseif n == 1
    L = [sym(1); 1 + alpha - x];
else
    L(1) = sym(1);
    L(2) = 1 + alpha - x;
    
    for i = 2:n
        Ai = - 1 / i;
        Bi = (2 * (i - 1) + alpha + 1) / i;
        Ci = (i + alpha - 1) / i;
        L(i + 1) = (Ai * x + Bi) * L(i) - Ci * L(i-1);
    end
    L = L';
end
end