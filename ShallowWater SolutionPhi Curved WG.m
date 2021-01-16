g = 980;
h0 = 2.5;
hg = 1.5;
delta = 2;
b = 60;
c = b + delta;
f = 5.0;
w = 2 * pi * f;

syms a k x r ;

J(a, k, x) = besselj(a, k * x);
DJ(a, k, x) = k * (J(a-1, k, x)-(a/k/x) *J(a, k, x));

Y(a, k, x) = bessely(a, k * x);
DY(a, k, x) = k * (Y(a-1, k, x)-(a/k/x) *Y(a, k, x));

RelDisp_0(k) = w^2 - g * k * h0 * k;
RelDisp_g(k) = w^2 - g * k * hg * k;

k0 = vpasolve(RelDisp_0 == 0, k, w/g/h0);
kg = vpasolve(RelDisp_g == 0, k, w/g/hg);

EqAlpha(a) = [J(a, k0, b) * DJ(a, kg, b) - DJ(a, k0, b) * J(a, kg, b)] * [Y(a, kg, c) * DY(a, k0, c) - Y(a, k0, c) * DY(a, kg, c)] - [DJ(a, k0, b) * Y(a, kg, b) - J(a, k0, b) * DY(a, kg, b)] * [DJ(a, kg, c) * Y(a, k0, c) - J(a, kg, c) * DY(a, k0, c)];

A = linspace(k0 * b, kg * b);

figure(1)
plot(A, EqAlpha(A))

alpha = vpasolve(EqAlpha == 0, a, [b * k0, b * kg]);

a = alpha(1)

Q = (J(a, kg, b) + Y(a, kg, b)) / J(a, k0, b);
P = (J(a, kg, c) + Y(a, kg, c)) / Y(a, k0, c);

sym r ;
phi(r) = piecewise(0<r<=b, Q * J(a, k0, r), b<r<c, J(a, kg, r) + Y(a, kg, r), r>=c, P * Y(a, k0, r));

r = b-5:0.01:c+5;

figure(2)
plot(r, phi(r), 'b-', 'LineWidth', 2)



