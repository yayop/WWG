%% ---------------- Propagation modes in curved WWG  ------------------- %%

%Author : Edgardo Rosas
%Project : Water Waveguiding (WWG)
%% ----------------- -------- Parameters ------------------------------- %%

% b : Inner curvature radius
% f : Frecuency
% omega : Angular freciency
% h0 : Fluid layer depth out waveguide
% hg : Fluid layer depth in waveguide
% g : Acceleration of gravity
% delta : Waveguide width

b = 5:2.5:100;
f = 2:0.5:5;
omega = 2*pi*f;
h0 = 2.5;
hg = 1.5;
g = 980;
delta = 10.0;
%% ---------------------- Waveguide band modes ------------------------- %%

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
%% ------------------------- Theoretical Data -------------------------- %%
%% ---------------------------- f = 2.0 hz ----------------------------- %%
R(1,:) = [3.54558,4.30395,5.06586,5.83011,6.59587,7.36255,8.12976,8.89724,9.66482,10.4324,11.1998,11.9671,12.7342,13.5011,14.2678,15.0343,15.8005,16.5665,17.3323,18.0979,18.8633,19.6285,20.3935,21.1583,21.923,22.6875,23.4518,24.216,24.9801,25.7441,26.5079,27.2716,28.0352,28.7988,29.5622,30.3256,31.0889,31.8521,32.6152];
I(1,:) = [0.867597,0.841382,0.814236,0.787256,0.760939,0.735493,0.710986,0.687417,0.664757,0.642964,0.621991,0.601792,0.582324,0.563544,0.545417,0.527906,0.510982,0.494615,0.478778,0.463449,0.448604,0.434224,0.42029,0.406783,0.393688,0.38099,0.368675,0.356729,0.345139,0.333894,0.322984,0.312396,0.302122,0.292152,0.282476,0.273087,0.263975,0.255133,0.246554];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(1)
title('$\Delta = 10.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(2)
title('$\Delta = 10.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 2.5 hz ----------------------------- %%
R(2,:) = [4.83149,5.81423,6.80126,7.7916,8.78446,9.7792,10.7753,11.7725,12.7705,13.769,14.768,15.7672,16.7667,17.7663,18.7661,19.7659,20.7658,21.7657,22.7657,23.7657,24.7656,25.7656,26.7657,27.7657,28.7657,29.7658,30.7659,31.7661,32.7662,33.7664,34.7667,35.767,36.7674,37.7678,38.7683,39.7688,40.7695,41.7702,42.771];
I(2,:) = [0.882681,0.847365,0.811014,0.774808,0.739415,0.705194,0.672316,0.640846,0.610784,0.582097,0.554736,0.528641,0.503752,0.480006,0.457346,0.435714,0.415058,0.395329,0.376481,0.358471,0.341258,0.324806,0.309079,0.294045,0.279673,0.265935,0.252802,0.240251,0.228256,0.216794,0.205844,0.195384,0.185396,0.17586,0.166758,0.158073,0.149788,0.141887,0.134355];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(4)
title('$\Delta = 10.0$ (cm), $f = 2.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 3.0 hz ----------------------------- %%
R(3,:) = [6.24651,7.47668,8.71165,9.9505,11.1925,12.4371,13.6838,14.9323,16.1822,17.4333,18.6854,19.9383,21.192,22.4462,23.7011,24.9563,26.2121,27.4682,28.7246,29.9814,31.2385,32.4959,33.7535,35.0115,36.2697,37.5282,38.7869,40.0459,41.3051,42.5646,43.8243,45.0843,46.3445,47.605,48.8657,50.1266,51.3878,52.6493,53.911];
I(3,:) = [0.927406,0.886284,0.844175,0.802169,0.760948,0.720934,0.682372,0.645394,0.610055,0.576362,0.544291,0.513799,0.484832,0.457329,0.431228,0.406465,0.382979,0.360709,0.339597,0.319588,0.300629,0.282669,0.265662,0.249561,0.234324,0.219908,0.206277,0.193391,0.181216,0.169719,0.158866,0.148626,0.138971,0.129872,0.121302,0.113234,0.105644,0.098508,0.0918029];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(5)
title('$\Delta = 10.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(6)
title('$\Delta = 10.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 3.5 hz ----------------------------- %%
R(4,:) = [7.82762,9.33493,10.8474,12.3642,13.8845,15.4079,16.9338,18.4619,19.9919,21.5236,23.0567,24.5912,26.1267,27.6634,29.2009,30.7393,32.2784,33.8183,35.3588,36.9,38.4417,39.984,41.5268,43.0701,44.6139,46.1581,47.7029,49.2481,50.7937,52.3397,53.8862,55.4331,56.9804,58.528,60.0761,61.6246,63.1734,64.7226,66.2722];
I(4,:) = [1.01978,0.976163,0.931335,0.886366,0.841935,0.798485,0.756302,0.715569,0.676395,0.638836,0.602914,0.568622,0.535935,0.504818,0.475226,0.44711,0.420417,0.395092,0.371081,0.348329,0.326784,0.306392,0.287104,0.26887,0.251642,0.235375,0.220023,0.205546,0.1919,0.179048,0.16695,0.15557,0.144873,0.134825,0.125393,0.116546,0.108253,0.100486,0.0932174];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(7)
title('$\Delta = 10.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(8)
title('$\Delta = 10.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 4.0 hz ----------------------------- %%
R(5,:) = [9.61679,11.4389,13.2666,15.0989,16.935,18.7743,20.6164,22.4609,24.3076,26.1562,28.0065,29.8584,31.7116,33.5661,35.4218,37.2785,39.1362,40.9949,42.8544,44.7146,46.5757,48.4374,50.2998,52.1629,54.0266,55.8909,57.7557,59.6211,61.487,63.3534,65.2204,67.0878,68.9557,70.8241,72.6929,74.5622,76.4319,78.3021,80.1727];
I(5,:) = [1.18369,1.14208,1.09844,1.05396,1.00938,0.96521,0.921778,0.879321,0.837998,0.797918,0.759152,0.721741,0.685706,0.651051,0.61777,0.585845,0.555253,0.525965,0.49795,0.471172,0.445595,0.421182,0.397895,0.375695,0.354544,0.334405,0.315241,0.297014,0.27969,0.263233,0.247608,0.232783,0.218725,0.205402,0.192784,0.18084,0.169542,0.158862,0.148771];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(9)
title('$\Delta = 10.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(10)
title('$\Delta = 10.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 4.5 hz ----------------------------- %%
R(6,:) = [11.6555,13.8382,16.0272,18.221,20.419,22.6203,24.8246,27.0314,29.2405,31.4516,33.6645,35.879,38.0949,40.3123,42.5308,44.7505,46.9713,49.193,51.4157,53.6392,55.8636,58.0886,60.3144,62.5409,64.768,66.9957,69.224,71.4528,73.6822,75.9121,78.1424,80.3733,82.6046,84.8363,87.0685,89.301,91.534,93.7674,96.0012];
I(6,:) = [1.44753,1.41477,1.37838,1.33982,1.30002,1.2596,1.21901,1.17856,1.13848,1.09893,1.06004,1.0219,0.98459,0.948153,0.912625,0.878033,0.84439,0.811705,0.77998,0.749211,0.719391,0.690511,0.662558,0.635517,0.609373,0.584107,0.559702,0.536138,0.513397,0.491459,0.470303,0.449909,0.430258,0.411329,0.393103,0.37556,0.35868,0.342445,0.326835];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(11)
title('$\Delta = 10.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(12)
title('$\Delta = 10.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ---------------------------- f = 5.0 hz ----------------------------- %%
R(7,:) = [13.9768,16.573,19.1761,21.7847,24.3979,27.0147,29.6348,32.2575,34.8827,37.51,40.1391,42.77,45.4024,48.0361,50.6712,53.3074,55.9447,58.583,61.2223,63.8624,66.5033,69.145,71.7874,74.4304,77.0741,79.7185,82.3633,85.0087,87.6547,90.3011,92.9479,95.5953,98.243,100.891,103.54,106.189,108.838,111.487,114.137];
I(7,:) = [1.83627,1.82192,1.80159,1.77713,1.74976,1.72033,1.68943,1.6575,1.62489,1.59183,1.55853,1.52515,1.4918,1.45858,1.42558,1.39287,1.36048,1.32847,1.29688,1.26573,1.23504,1.20483,1.17513,1.14593,1.11725,1.08909,1.06146,1.03436,1.00779,0.981745,0.956231,0.931242,0.906773,0.882821,0.85938,0.836446,0.814011,0.792071,0.770617];
% ---------------------- Figure: alpha_R, alpha_I ---------------------- %%
figure(13)
title('$\Delta = 10.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
% ------------------ Figure: alpha_R / b, alpha_I / b ------------------ %%
figure(14)
title('$\Delta = 10.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% ----------------------------- 3D PLOT ------------------------------- %%
bMat = repmat(b', 1, numel(f));
fMat = repmat(f, numel(b), 1);
figure(16)
plot3(bMat, fMat, R./bMat', 'b*-', 'LineWidth', 2)
grid on
zlabel('Re($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 10.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(17)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 10.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(18)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 2.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
set(gca,'zscale','log')
set(gca,'xscale','log')