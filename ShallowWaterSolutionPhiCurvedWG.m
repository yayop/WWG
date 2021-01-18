g = 980.0;
h0 = 2.5;
hg = 1.5;
delta = 10.0;
b = 800.0;
c = b + delta;
f = 5.0;
w = 2 * pi * f;

ri = 750;
rf = 850;

syms a k x r ;

J(a, k, x) = 10*besselj(a, k * x);
DJ(a, k, x) = k * (J(a-1, k, x)-(a/k/x) *J(a, k, x));

Y(a, k, x) = 10*bessely(a, k * x);
DY(a, k, x) = k * (Y(a-1, k, x)-(a/k/x) *Y(a, k, x));

RelDisp_0(k) = w^2 - g * k * tanh(h0 * k);
RelDisp_g(k) = w^2 - g * k * tanh(hg * k);

k0 = vpasolve(RelDisp_0 == 0, k, w/sqrt(g*h0));
k0 = k0(1);

kg = vpasolve(RelDisp_g == 0, k,  w/sqrt(g*hg));
kg = kg(1);

MAT(a) = [J(a, k0, b) -J(a, kg, b) -Y(a, kg, b) 0; 0, J(a, kg, c) Y(a, kg, c) -Y(a, k0, c); b*DJ(a, k0, b) -b*DJ(a, kg, b) -b*DY(a, kg, b) 0; 0 c*DJ(a, kg, c) c*DY(a, kg, c) -c*DY(a, k0, c)]
EqAlpha(a) = det(MAT(a))

A = linspace(k0 * b, kg * b);

figure(1)
plot(A, EqAlpha(A))

alpha = vpasolve(EqAlpha == 0, a, [b * k0, b * kg]);

a = alpha(1)

Q = (J(a, kg, b) + Y(a, kg, b)) / J(a, k0, b);
P = (J(a, kg, c) + Y(a, kg, c)) / Y(a, k0, c);

sym r ;

phi(r) = piecewise(0<r<=b, Q * J(a, k0, r), b<r<c, +J(a, kg, r) + Y(a, kg, r), r>=c, P * Y(a, k0, r));

r = linspace(ri,rf,500);

figure(2)
plot(r, phi(r), 'b-', 'LineWidth', 2)
grid on
xlabel('$r$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('$\eta$ (cm)', 'interpreter', 'latex','Fontsize',14)


figure(3)
r = linspace(ri,rf);
psi= linspace(0,pi/2);
[R PSI] = meshgrid(r, psi);
Z = phi(R).*cos(a*PSI);
Z = double(Z);
surf(R.*cos(PSI), R.*sin(PSI), Z)
c = colorbar('LineWidth',1.5,'FontSize',13);
colormap(jet)

k0<a/b<kg