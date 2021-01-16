%% ---------------- Propagation modes in curved WWG  ------------------- %%
%Author : Edgardo Rosas
%Project : Water Waveguiding (WWG)
%% Parameters
% b : Inner curvature radii
% f : Frecuency
% omega : Angular freciency
% h0 : Fluid layer depth out waveguide
% hg : Fluid layer depth in waveguide

b = 5:2.5:100;
f = 2:0.5:5;
omega = 2*pi*f;
g = 980;
h0 = 2.5;
hg = 1.5;
%% Waveguide band modes
syms k w
RelDisp_0(w, k) = w^2-g*k*tanh(k*h0);
RelDisp_g(w, k) = w^2-g*k*tanh(k*hg);

K0 = [];
Kg = [];
for i = 1:1:numel(f)
    uk0 = vpasolve(RelDisp_0(omega(i), k) == 0, k, omega(i)/sqrt(g*h0));
    K0(i) = uk0(1);
    
    ukg = vpasolve(RelDisp_g(omega(i), k) == 0, k, omega(i)/sqrt(g*hg));
    Kg(i) = ukg(1);
    
    kg(i,:) = repmat(Kg(i),1,numel(b));
    k0(i,:) = repmat(K0(i),1,numel(b));
end
%% Data
%% f = 2.0 hz, Delta = 2
R(1,:) = [1.18202,1.87945,2.57516,3.27044,3.96562,4.66082,5.35604,6.05129,6.74656,7.44183,8.13708,8.83232,9.52753,10.2227,10.9179,11.6129,12.308,13.003,13.6979,14.3928,15.0877,15.7825,16.4772,17.1719,17.8665,18.5611,19.2557,19.9501,20.6446,21.3389,22.0332,22.7275,23.4217,24.1159,24.81,25.5041,26.1981,26.8921,27.5861];
I(1,:) = [1.11156,1.20974,1.28485,1.34556,1.39638,1.44,1.47812,1.5119,1.54218,1.56957,1.59451,1.61738,1.63846,1.65797,1.67611,1.69304,1.70887,1.72373,1.73771,1.75089,1.76334,1.77512,1.78629,1.7969,1.80698,1.81658,1.82572,1.83445,1.84279,1.85076,1.85838,1.86568,1.87268,1.87939,1.88583,1.89201,1.89795,1.90365,1.90914];
% Figure: alpha_R, alpha_I
figure(1)
title('$\Delta = 2.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(1,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(1,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(2)
title('$\Delta = 2.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(1,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(1,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% f = 2.5 hz, Delta = 1
R(2,:) = [1.75863,2.67603,3.59124,4.50551,5.41923,6.33256,7.24555,8.15824,9.07066,9.98282,10.8947,11.8064,12.7179,13.6291,14.5402,15.4511,16.3617,17.2722,18.1825,19.0927,20.0027,20.9125,21.8222,22.7317,23.6412,24.5504,25.4596,26.3686,27.2775,28.1862,29.0949,30.0034,30.9119,31.8202,32.7284,33.6366,34.5446,35.4525,36.3604];
I(2,:) = [1.13912,1.22638,1.29197,1.34408,1.38699,1.42322,1.45439,1.48158,1.50556,1.5269,1.54603,1.56329,1.57893,1.59317,1.60618,1.61811,1.62908,1.63919,1.64851,1.65714,1.66513,1.67255,1.67943,1.68582,1.69176,1.69729,1.70243,1.70722,1.71168,1.71582,1.71968,1.72326,1.72659,1.72968,1.73255,1.73521,1.73766,1.73993,1.74202];
% Figure: alpha_R, alpha_I
figure(3)
title('$\Delta = 1.0$ (cm), $f = 2.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(2,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(2,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(4)
title('$\Delta = 2.0$ (cm), $f = 2.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(2,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(2,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% f = 3.0 hz, Delta = 2
R(3,:) = [2.3931,3.55794,4.72031,5.88137,7.04153,8.20095,9.35974,10.518,11.6757,12.8329,13.9897,15.1461,16.3021,17.4577,18.613,19.768,20.9226,22.077,23.2311,24.3849,25.5385,26.6918,27.8449,28.9978,30.1505,31.303,32.4553,33.6073,34.7593,35.911,37.0625,38.2139,39.3652,40.5163,41.6672,42.818,43.9686,45.1192,46.2695];
I(3,:) = [1.18044,1.2583,1.31584,1.36076,1.39709,1.42719,1.45259,1.47431,1.49306,1.50939,1.5237,1.53629,1.54742,1.55728,1.56604,1.57383,1.58075,1.58691,1.59238,1.59724,1.60154,1.60532,1.60865,1.61156,1.61407,1.61624,1.61807,1.6196,1.62085,1.62184,1.62259,1.62311,1.62342,1.62354,1.62347,1.62323,1.62283,1.62227,1.62158];
% Figure: alpha_R, alpha_I
figure(5)
title('$\Delta = 2.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(3,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(3,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(6)
title('$\Delta = 2.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(3,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(3,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% f = 3.5 hz, Delta = 2
R(4,:) = [3.10131,4.5483,5.99312,7.43661,8.87907,10.3206,11.7614,13.2015,14.641,16.0798,17.518,18.9558,20.393,21.8298,23.2662,24.7021,26.1377,27.5729,29.0078,30.4423,31.8765,33.3104,34.7441,36.1774,37.6105,39.0434,40.476,41.9083,43.3405,44.7724,46.2041,47.6356,49.067,50.4981,51.9291,53.3599,54.7905,56.2209,57.6512];
I(4,:) = [1.25033,1.32192,1.37428,1.41461,1.44674,1.47294,1.49466,1.51286,1.52826,1.54136,1.55254,1.56211,1.57031,1.57732,1.58329,1.58836,1.59263,1.5962,1.59912,1.60149,1.60334,1.60473,1.6057,1.60629,1.60653,1.60645,1.60608,1.60544,1.60456,1.60344,1.60211,1.60059,1.59888,1.59701,1.59498,1.5928,1.59048,1.58804,1.58548];
% Figure: alpha_R, alpha_I
figure(7)
title('$\Delta = 2.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(4,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(4,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(8)
title('$\Delta = 2.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(4,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(4,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% f = 4.0 hz, Delta = 5
R(5,:) = [3.89862,5.67007,7.44061,9.2104,10.9795,12.7478,14.5155,16.2825,18.0489,19.8147,21.58,23.3447,25.1089,26.8727,28.636,30.3989,32.1614,33.9235,35.6852,37.4466,39.2077,40.9684,42.7289,44.489,46.2489,48.0085,49.7678,51.5269,53.2858,55.0444,56.8028,58.561,60.319,62.0768,63.8344,65.5918,67.349,69.1061,70.8629];
I(5,:) = [1.36666,1.43707,1.48872,1.5285,1.56014,1.58585,1.60708,1.62478,1.63966,1.65223,1.66288,1.67189,1.67952,1.68595,1.69134,1.69582,1.6995,1.70245,1.70477,1.70652,1.70776,1.70853,1.70888,1.70885,1.70847,1.70776,1.70677,1.70551,1.70399,1.70225,1.7003,1.69816,1.69583,1.69333,1.69068,1.68788,1.68495,1.68189,1.67872];
% Figure: alpha_R, alpha_I
figure(9)
title('$\Delta = 2.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(5,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(5,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(10)
title('$\Delta = 2.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(5,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(5,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% f = 4.5 hz, Delta = 2
R(6,:) = [4.7973,6.94164,9.08785,11.2348,13.382,15.5291,17.676,19.8225,21.9688,24.1147,26.2602,28.4053,30.55,32.6944,34.8385,36.9822,39.1255,41.2686,43.4113,45.5538,47.6959,49.8378,51.9794,54.1208,56.2619,58.4028,60.5434,62.6838,64.824,66.964,69.1038,71.2434,73.3828,75.522,77.6611,79.7999,81.9386,84.0772,86.2155];
I(6,:) = [1.54989,1.62634,1.68356,1.72843,1.76472,1.79471,1.81989,1.84127,1.85959,1.87539,1.88907,1.90097,1.91132,1.92034,1.9282,1.93504,1.94097,1.94609,1.95048,1.95423,1.95738,1.96,1.96213,1.96382,1.9651,1.96601,1.96657,1.96681,1.96676,1.96644,1.96587,1.96506,1.96403,1.96279,1.96137,1.95976,1.95799,1.95607,1.95399];
% Figure: alpha_R, alpha_I
figure(11)
title('$\Delta = 2.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(6,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(6,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(12)
title('$\Delta = 2.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(6,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(6,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';

%% f = 5.0 hz, Delta = 2
R(7,:) = [5.80469,8.37306,10.9482,13.5268,16.1075,18.6894,21.272,23.855,26.4383,29.0217,31.6051,34.1885,36.7719,39.3551,41.9382,44.5212,47.1041,49.6868,52.2694,54.8518,57.434,60.0161,62.598,65.1798,67.7614,70.3428,72.9241,75.5053,78.0863,80.6672,83.2479,85.8285,88.409,90.9893,93.5695,96.1495,98.7295,101.309,103.889];
I(7,:) = [1.81927,1.9109,1.98153,2.03847,2.08579,2.12596,2.16061,2.19087,2.21755,2.24126,2.26247,2.28153,2.29874,2.31434,2.32852,2.34144,2.35324,2.36403,2.37391,2.38297,2.39128,2.39891,2.40592,2.41234,2.41824,2.42364,2.42859,2.43312,2.43726,2.44102,2.44445,2.44755,2.45035,2.45287,2.45512,2.45711,2.45888,2.46042,2.46174];
figure(13)
title('$\Delta = 2.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(7,:), 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(7,:), 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%
figure(14)
title('$\Delta = 2.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
yyaxis left 
plot(b, R(7,:)./b, 'b*-', 'LineWidth', 2)
grid on
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
ylabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
yyaxis right
plot(b, I(7,:)./b, 'r*-', 'LineWidth', 2)
ylabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
%% 3D
bMat = repmat(b', 1, numel(f));
fMat = repmat(f, numel(b), 1);
figure(16)
plot3(bMat, fMat, R./bMat', 'b*-', 'LineWidth', 2)
grid on
zlabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 2.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(17)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 2.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(18)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 2.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
set(gca,'zscale','log')
set(gca,'xscale','log')
