i_1 = 16;
i_2 = 31;

syms x y b w kg k0
l_g(kg, b) = sqrt(kg^2-b^2);
l_0(k0, b) = sqrt(b^2-k0^2);

A(k0, kg, b) = cos(l_g(kg, b) * a/2) * exp(l_0(k0, b)* a/2);
B(k0, kg, b) = -sin(l_g(kg, b) * a/2) * exp(l_0(k0, b)* a/2);

phi_s(x, k0, kg, b) = piecewise(abs(x) >= a/2, A(k0, kg, b) * exp(- l_0(k0, b) * abs(x)), abs(x) < a/2, cos(l_g(kg, b) * x));
phi_a(x, k0, kg, b) = piecewise(x < -a/2, B(k0, kg, b) * exp(l_0(k0, b) * x), abs(x) < a/2, sin(l_g(kg, b) * x), x > a/2, - B(k0, kg, b) * exp(-l_0(k0, b) * x));
I_a(k0, kg, b) = -int(phi_a(x,k0, kg, b), x, -10*a,0)+int(+phi_a(x,k0, kg, b), x, 0,10*a);
I_s(k0, kg, b) = int(phi_s(x, k0, kg, b), x, -10*a,10*a);
PHI(x, y) = phi_s(x, K0(i_1), Kg(i_1), BetaFWs(i_1)) * cos(BetaFWs(i_1)*y);
figure(1)
%fsurf(PHI, [-4*a, 4*a, 0, 30])
figure(2)
x = -2.5*a:0.05:+2.5*a;
subplot(2,1,1);
p16_s = I_s(K0(i_2), Kg(i_2), BetaFWs(i_2))/I_s(K0(i_1), Kg(i_1), BetaFWs(i_1));
plot(x/a, p16_s * phi_s(x,K0(i_1), Kg(i_1), BetaFWs(i_1)), 'b-', 'LineWidth', 2)
grid on
hold on
plot(x/a, phi_s(x,K0(i_2), Kg(i_2), BetaFWs(i_2)), 'r-.', 'LineWidth', 2)
hold off
xlabel('$x/a$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\eta/\eta^*$', 'interpreter', 'latex', 'fontsize', 14)
l = legend('$f = 3.5$ (hz)', '$f=5.0$ (hz)');
set(l, 'interpreter', 'latex')

subplot(2,1,2);
p16_a = I_a(K0(i_2), Kg(i_2), BetaFWa(i_2))/I_a(K0(i_1), Kg(i_1), BetaFWa(i_1));
plot(x/a, p16_a * phi_a(x,K0(i_1), Kg(i_1), BetaFWa(i_1)),'b-', 'LineWidth', 2)
grid on
hold on
plot(x/a, phi_a(x,K0(i_2), Kg(i_2), BetaFWa(i_2)),'r-.', 'LineWidth', 2)
%plot(x/a, sin(l_g(omega(31), betaSWa(31)) * x),'g-', 'LineWidth', 3)
hold off
xlabel('$x/a$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\eta/\eta^*$', 'interpreter', 'latex', 'fontsize', 14)
l = legend('$f = 3.5$ (hz)', '$f=5.0$ (hz)');
set(l,'interpreter', 'latex')
