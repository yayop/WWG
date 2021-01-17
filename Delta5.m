%% ---------------- Propagation modes in curved WWG  ------------------- %%
%Author : Edgardo Rosas
%Project : Water Waveguiding (WWG)
%% Parameters
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
delta = 5.0;
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
%% f = 2.0 hz, Delta = 5
R(1,:) = [2.06287,2.80903,3.55221,4.29314,5.03231,5.77002,6.5065,7.24191,7.97637,8.70999,9.44284,10.175,10.9065,11.6375,12.3678,13.0977,13.8271,14.5561,15.2847,16.0128,16.7406,17.4681,18.1952,18.922,19.6485,20.3748,21.1007,21.8264,22.5519,23.2771,24.0021,24.7268,25.4514,26.1758,26.8999,27.6239,28.3477,29.0713,29.7948];
I(1,:) = [0.930588,0.945404,0.955795,0.96275,0.967098,0.969445,0.970233,0.969784,0.968339,0.966084,0.963162,0.959685,0.955744,0.951413,0.946751,0.941809,0.936627,0.931242,0.925682,0.919974,0.914138,0.908194,0.902158,0.896043,0.889864,0.88363,0.877351,0.871035,0.864691,0.858324,0.851941,0.845547,0.839147,0.832746,0.826346,0.819952,0.813566,0.807191,0.800831];
% Figure: alpha_R, alpha_I
figure(1)
title('$\Delta = 5.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 2.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% f = 2.5 hz, Delta = 5
R(2,:) = [2.89507,3.86618,4.83656,5.80558,6.77317,7.73939,8.70433,9.66811,10.6308,11.5926,12.5534,13.5134,14.4726,15.4312,16.389,17.3463,18.303,19.2591,20.2148,21.1699,22.1246,23.0789,24.0328,24.9864,25.9395,26.8924,27.8449,28.7971,29.749,30.7006,31.652,32.6031,33.554,34.5046,35.4551,36.4053,37.3553,38.3051,39.2547];
I(2,:) = [0.953625,0.94946,0.943297,0.935643,0.926892,0.917338,0.907194,0.896622,0.88574,0.874639,0.863389,0.852043,0.840644,0.829226,0.817815,0.806434,0.7951,0.783827,0.772626,0.761508,0.750479,0.739546,0.728715,0.717988,0.70737,0.696863,0.686469,0.67619,0.666026,0.65598,0.64605,0.636238,0.626544,0.616967,0.607508,0.598165,0.588939,0.579829,0.570833];
% Figure: alpha_R, alpha_I
figure(3)
title('$\Delta = 5.0$ (cm), $f = 2.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 2.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% f = 3.0 hz, Delta = 5
R(3,:) = [3.81262,5.03121,6.25215,7.47335,8.69402,9.91383,11.1327,12.3505,13.5674,14.7833,15.9983,17.2125,18.426,19.6387,20.8507,22.0621,23.2729,24.4832,25.6929,26.9021,28.1109,29.3192,30.5271,31.7346,32.9417,34.1485,35.355,36.5611,37.767,38.9726,40.1779,41.3829,42.5878,43.7924,44.9967,46.2009,47.4049,48.6087,49.8123];
I(3,:) = [1.00383,0.985216,0.965725,0.946084,0.926541,0.907204,0.88813,0.869352,0.850887,0.832746,0.814933,0.797451,0.780297,0.763471,0.746969,0.730786,0.714917,0.699358,0.684103,0.669147,0.654484,0.640109,0.626015,0.612198,0.598653,0.585374,0.572356,0.559593,0.547082,0.534817,0.522793,0.511006,0.499452,0.488126,0.477023,0.46614,0.455473,0.445018,0.43477];
% Figure: alpha_R, alpha_I
figure(5)
title('$\Delta = 5.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 3.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% f = 3.5 hz, Delta = 5
R(4,:) = [4.83845,6.33334,7.83406,9.33735,10.8416,12.346,13.8501,15.3537,16.8567,18.359,19.8606,21.3615,22.8618,24.3614,25.8604,27.3589,28.8568,30.3542,31.851,33.3475,34.8434,36.339,37.8342,39.329,40.8235,42.3176,43.8114,45.3049,46.7982,48.2912,49.7839,51.2764,52.7687,54.2608,55.7527,57.2444,58.7359,60.2273,61.7185];
I(4,:) = [1.09661,1.06963,1.04117,1.01291,0.985335,0.958599,0.932716,0.90766,0.88339,0.859864,0.837043,0.814888,0.793367,0.772448,0.752104,0.732309,0.713041,0.694279,0.676004,0.658197,0.640843,0.623926,0.607432,0.591347,0.575659,0.560357,0.545429,0.530865,0.516654,0.502787,0.489256,0.476051,0.463164,0.450588,0.438315,0.426337,0.414648,0.403241,0.392109];
% Figure: alpha_R, alpha_I
figure(7)
title('$\Delta = 5.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 3.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
R(5,:) = [5.99832,7.80602,9.62282,11.4449,13.27,15.0967,16.9242,18.7521,20.58,22.4077,24.2352,26.0622,27.8889,29.7152,31.541,33.3664,35.1913,37.0159,38.84,40.6638,42.4871,44.3102,46.1328,47.9552,49.7772,51.5989,53.4204,55.2416,57.0625,58.8832,60.7036,62.5238,64.3439,66.1637,67.9833,69.8028,71.622,73.4412,75.2601];
I(5,:) = [1.25148,1.22442,1.19331,1.16126,1.1295,1.09851,1.06847,1.03941,1.01133,0.984169,0.957894,0.932454,0.907803,0.883897,0.860698,0.83817,0.81628,0.794997,0.774297,0.754152,0.734541,0.715443,0.696839,0.67871,0.66104,0.643813,0.627015,0.610633,0.594653,0.579064,0.563854,0.549012,0.534529,0.520394,0.506599,0.493135,0.479992,0.467164,0.454642];
% Figure: alpha_R, alpha_I
figure(9)
title('$\Delta = 5.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 4.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
%% f = 4.5 hz, Delta = 5
R(6,:) = [7.31736,9.48271,11.6596,13.8444,16.0343,18.2279,20.4238,22.6212,24.8197,27.0188,29.2183,31.418,33.6177,35.8174,38.017,40.2165,42.4157,44.6147,46.8135,49.012,51.2103,53.4083,55.606,57.8036,60.0008,62.1979,64.3946,66.5912,68.7875,70.9836,73.1796,75.3753,77.5708,79.7661,81.9612,84.1561,86.3509,88.5455,90.74];
I(6,:) = [1.49072,1.47465,1.45023,1.42202,1.39229,1.36219,1.3323,1.30294,1.27422,1.24623,1.21897,1.19244,1.16661,1.14147,1.11699,1.09313,1.06988,1.0472,1.02508,1.00348,0.982392,0.961789,0.941655,0.921972,0.902724,0.883895,0.865473,0.847443,0.829793,0.812512,0.795588,0.779012,0.762773,0.746863,0.731272,0.715993,0.701017,0.686337,0.671945];
% Figure: alpha_R, alpha_I
figure(11)
title('$\Delta = 5.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 4.5$ (hz)', 'interpreter', 'latex','Fontsize',14)
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

%% f = 5.0 hz, Delta = 5
R(7,:) = [8.81522,11.3903,13.9787,16.5767,19.1819,21.7926,24.4073,27.0251,29.6452,32.2671,34.8904,37.5147,40.1398,42.7654,45.3915,48.0179,50.6445,53.2712,55.898,58.5248,61.1517,63.7785,66.4052,69.0319,71.6585,74.2849,76.9113,79.5375,82.1636,84.7896,87.4154,90.0411,92.6666,95.292,97.9173,100.542,103.167,105.792,108.417];
I(7,:) = [1.83303,1.84138,1.83617,1.82297,1.80503,1.7843,1.76198,1.7388,1.71524,1.69158,1.66801,1.64465,1.62157,1.59882,1.57642,1.55439,1.53272,1.51143,1.4905,1.46993,1.44971,1.42983,1.41029,1.39107,1.37217,1.35358,1.33528,1.31727,1.29954,1.28209,1.2649,1.24798,1.2313,1.21487,1.19868,1.18273,1.167,1.1515,1.13621];
figure(13)
title('$\Delta = 5.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm), $f = 5.0$ (hz)', 'interpreter', 'latex','Fontsize',14)
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
title('$\Delta = 5.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(17)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 5.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
figure(18)
plot3(bMat, fMat, I./bMat', 'r*-', 'LineWidth', 2)
grid on
zlabel('Im($\alpha/b$)', 'interpreter', 'latex','Fontsize',16)
ylabel('$f$ (hz)', 'interpreter', 'latex','Fontsize',14)
xlabel('$b$ (cm)', 'interpreter', 'latex','Fontsize',14)
title('$\Delta = 2.0$ (cm)', 'interpreter', 'latex','Fontsize',14)
set(gca,'zscale','log')
set(gca,'xscale','log')