i_1 = 16;
i_2 = 31;

syms x y b w
l_g(w, b) = sqrt(w^2/(g*hg)-b^2);
l_0(w, b) = sqrt(b^2-w^2/(g*h0));

A(w, b) = cos(l_g(w, b) * a/2) * exp(l_0(w, b)* a/2);
B(w, b) = -sin(l_g(w, b) * a/2) * exp(l_0(w, b)* a/2);

phi_s(x, w, b) = piecewise(abs(x) >= a/2, A(w, b) * exp(- l_0(w, b) * abs(x)), abs(x) < a/2, cos(l_g(w, b) * x));
phi_a(x, w, b) = piecewise(x < -a/2, B(w, b) * exp(l_0(w, b) * x), abs(x) < a/2, sin(l_g(w, b) * x), x > a/2, - B(w, b) * exp(-l_0(w, b) * x));
I_a(w,b) = -int(phi_a(x,w,b), x, -10*a,0)+int(+phi_a(x,w,b), x, 0,10*a);
I_s(w,b) = int(phi_s(x,w,b), x, -10*a,10*a);

x = -2.5*a:0.05:+2.5*a;

figure(1)
subplot(2,1,1);
p16_s = I_s(omega(i_2),betaSWs(i_2))/I_s(omega(i_1),betaSWs(i_1));
plot(x/a, p16_s * phi_s(x,omega(i_1),betaSWs(i_1)), 'b-', 'LineWidth', 2)
grid on
hold on
plot(x/a, phi_s(x,omega(i_2),betaSWs(i_2)), 'r-.', 'LineWidth', 2)
hold off
xlabel('$x/a$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\eta/\eta^*$', 'interpreter', 'latex', 'fontsize', 14)
l = legend('$f = 3.5$ (hz)', '$f=5.0$ (hz)');
set(l, 'interpreter', 'latex')

subplot(2,1,2);
p16_a = I_a(omega(i_2),betaSWa(i_2))/I_a(omega(i_1),betaSWa(i_1));
plot(x/a, p16_a * phi_a(x,omega(i_1),betaSWa(i_1)),'b-', 'LineWidth', 2)
grid on
hold on
plot(x/a, phi_a(x,omega(i_2),betaSWa(i_2)),'r-.', 'LineWidth', 2)
%plot(x/a, sin(l_g(omega(31), betaSWa(31)) * x),'g-', 'LineWidth', 3)
hold off
xlabel('$x/a$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\eta/\eta^*$', 'interpreter', 'latex', 'fontsize', 14)
l = legend('$f = 3.5$ (hz)', '$f=5.0$ (hz)');
set(l,'interpreter', 'latex')

