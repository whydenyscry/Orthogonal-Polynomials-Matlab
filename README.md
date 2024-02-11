# Classical Orthogonal Polynomials
## General Orthogonal Polynomials. Definitions
**Orthogonality on intervals.** A set of polynomials $\lbrace p_n(x)\rbrace_{n=0}^{\infty}$ is said to be orthogonal on $\left(a,b\right)$ with respect to the weight function $\omega\left(x\right)\geq0$ if
$$\int_a^bp_n\left(x\right)p_m\left(x\right)\omega\left(x\right)\mathrm{d}x=0,\quad n\ne m.$$

**Orthonormality on intervals.** A set of polynomials $\left\lbrace p_n\left(x\right)\right\rbrace_{n=0}^{\infty}$ is said to be orthonormal on $\left(a,b\right)$ with respect to the weight function $\omega\left(x\right)\geq0$ if
$$\int_a^bp_n\left(x\right)p_m\left(x\right)\omega\left(x\right)\mathrm{d}x=\delta_{nm},$$
where $\delta_{nm}$ is Kronecker delta.

**Recurrence relations.** Assume that $p_{-1}\left(x\right)\equiv0$, then
\begin{equation}
	p_{n+1}\left(x\right)=\left(A_nx+B_n\right)p_n\left(x\right)-C_np_{n-1}\left(x\right),
\end{equation}
here $A_n,B_n\left(n\geq 0\right)$, and $C_n\left(n\geq 1\right)$ are real constants.

**Rodrigues' formula.** Orthogonal polynomials can be expressed through Rodrigue's formula, which gives an analytic expression for polynomials through derivatives:
$$p_{n}(x)=\frac{1}{\kappa_{n}\omega(x)}\frac{{\mathrm{d}}^{n}}{{\mathrm{d}x}^{n}}\left[\omega(x)(F(x))^{n}\right].$$

## Jacobi polynomials
The Jacobi polynomials $p_n\left(x\right)=P_{n}^{(\alpha ,\beta )}\left(x\right)$ are a class of orthogonal polynomials orthogonal on an interval $\left(-1,1\right)$ with a weight function $\omega\left(x\right)=\left(1-x\right)^\alpha\left(1+x\right)^\beta$. Gegenbauer, Chebyshev polynomials of all kinds and Legendre polynomials are special cases of Jacobi polynomials.

**Definition.** For $z\in\mathbb{C}$ Jacobi polynomials can be defined as
$$P_{n}^{(\alpha ,\beta )}(z)={\frac {\Gamma (\alpha +n+1)}{n!\,\Gamma (\alpha +\beta +n+1)}}\sum _{m=0}^{n}{n \choose m}{\frac {\Gamma (\alpha +\beta +n+m+1)}{\Gamma (\alpha +m+1)}}\left({\frac {z-1}{2}}\right)^{m}.$$

For $x\in\mathbb{R}$ Jacobi polynomials can be defined as
$$P_{n}^{(\alpha ,\beta )}(x)=\sum _{s=0}^{n}{n+\alpha  \choose n-s}{n+\beta  \choose s}\left({\frac {x-1}{2}}\right)^{s}\left({\frac {x+1}{2}}\right)^{n-s}.$$

Another representation can be obtained using the Rodrigues' formula:
$$P_{n}^{(\alpha ,\beta )}(x)=\frac{1}{\left(-2\right)^nn!}\left(1-x\right)^{-\alpha}\left(1+x\right)^{-\beta}\frac{{\mathrm{d}}^{n}}{{\mathrm{d}x}^{n}}\left[\left(1-x\right)^\alpha\left(1+x\right)^\beta\left(1-x^2\right)^{n}\right],$$
here for Jacobi polynomials $\kappa_{n}=\left(-2\right)^nn!, F\left(x\right)=\left(1-x^2\right).$

**Recurrence relations.** Coefficients are defined as:
\begin{equation}
P_{n+1}^{(\alpha ,\beta )}(x)=(A_{n}x+B_{n})P_{n}^{(\alpha ,\beta )}-C_{n}P_{n-1}^{(\alpha ,\beta )},
\end{equation}
where 

