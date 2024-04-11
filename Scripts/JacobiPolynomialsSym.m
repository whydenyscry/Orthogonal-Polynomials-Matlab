function P = JacobiPolynomialsSym(alpha,beta,n,x)

A0 = 1/2*(alpha+beta)+1;
B0 = 1/2*(alpha-beta);

P(1) = sym(1);
P(2) = A0*x+B0;

for i = 2:n
       Ai = (2*(i-1)+alpha+beta+1)*(2*(i-1)+alpha+beta+2)/2/((i-1)+1)/((i-1)+alpha+beta+1);
       Bi = (alpha^2-beta^2)*(2*(i-1)+alpha+beta+1)/2/((i-1)+1)/((i-1)+alpha+beta+1)/(2*(i-1)+alpha+beta);
       Ci = ((i-1)+alpha)*((i-1)+beta)*(2*(i-1)+alpha+beta+2)/((i-1)+1)/((i-1)+alpha+beta+1)/(2*(i-1)+alpha+beta);

       P(i+1)=(Ai*x+Bi)*P(i)-Ci*P(i-1);
end
end