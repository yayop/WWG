g = 980;
h0 = 2.5;
hg = 1.5;
f = 2:0.1:5;
omega = 2*pi*f;
a = 2.5;

syms k w h b
RELDISP(k, w, h) = -g*k.*tanh(k*h) + w^2;
RELDISPC(k, w, h) = g*k.*tan(k*h) + w^2;

%SHALLOW WATER
kg(w) = w/sqrt(g*hg);
k0(w) = w/sqrt(g*h0);

SWs(w, b) = sqrt(kg(w)^2-b^2)*tan((a/2)*sqrt(kg(w)^2-b^2))-sqrt(b^2-k0(w)^2);
SWa(w, b) = sqrt(kg(w)^2-b^2)*cot((a/2)*sqrt(kg(w)^2-b^2))+sqrt(b^2-k0(w)^2);

betaSWs = [];
betaSWa = [];
K0 = [];
Kg = [];
N_w = length(omega);
N_j = 5;
KK = [];
K = [];
betaFWs = [];
betaFWa = [];
BetaFWs = [];
BetaFWa = [];
beta_new = [];
for i = 1:1:length(omega)
    %BETA SHALLOW WATER SOLUCION SIMETRICA
    SWfs(b) = SWs(omega(i),b);
    bbs = vpasolve(SWfs == 0, b, [k0(omega(i))+eps, kg(omega(i))-eps]);
    if isempty(bbs)
        betaSWs(i) = NaN;
    else
        betaSWs(i) = bbs(1);
    end
    %BETA SHALLOW WATER SOLUCION ANTI-SIMETRICA
    SWfa(b) = SWa(omega(i),b);
    bba = vpasolve(SWfa == 0, b, [k0(omega(i))+eps, kg(omega(i))-eps]);
    if isempty(bba)
        betaSWa(i) = NaN;
    else
        betaSWa(i) = bba(1);
    end
    %CALCULO VECTOR DE ONDA Z SOL IMAGINARIA FUERA DE LA GUIA
    RELDISP0f(k) =  RELDISP(k, omega(i), h0);
    uk0 = vpasolve(RELDISP0f == 0, k, 1);
    K0(i) = uk0(1);
    %CALCULO VECTOR DE ONDA Z SOL IMAGINARIA Z EN LA GUIA
    RELDISPgf(k) = RELDISP(k, omega(i), hg);
    ukg = vpasolve(RELDISPgf == 0, k, 1);
    Kg(i) = ukg(1);
    %SIMETRICA
    W(b) = sqrt(ukg^2-b^2).*tan((a/2)*sqrt(ukg^2-b^2))-sqrt(b^2-uk0^2);
    Bfws = vpasolve(W == 0, b, [uk0+eps, ukg-eps]);
    if isempty(Bfws)
        BetaFWs(i) = NaN;
    else
        BetaFWs(i) = Bfws(1);
    end
    %ANTISIMETRICA
    Y(b) = sqrt(ukg^2-b^2).*cot((a/2)*sqrt(ukg^2-b^2))+sqrt(b^2-uk0^2);
    Bfwa = vpasolve(Y == 0, b, [uk0+eps, ukg-eps]);
    if isempty(Bfwa)
        BetaFWa(i) = NaN;
    else
        BetaFWa(i) = Bfwa(1);
    end
    for j = 1:2:N_j
        %CALCULO VECTOR DE ONDA Z SOL REAL FUERA DE LA GUIA
        RELDISPC0f(k) =  RELDISPC(k, omega(i),h0);
        uk = vpasolve(RELDISPC0f == 0, k, [j*pi/2/h0+0.1, (j+2)*pi/2/h0-0.1]);
        K(i,j) = uk(1);
        %CALCULO VECTOR DE ONDA Z SOL REAL EN LA GUIA
        RELDISPCgf(k) = RELDISPC(k, omega(i),hg);
        ukk = vpasolve(RELDISPCgf == 0, k, [j*pi/2/hg+0.1, (j+2)*pi/2/hg-0.1]);
        KK(i,j) = ukk(1);
        %TERMINOS EQ TRASENDENTAL PARA BETA SIMETRICO
        R(b) = sqrt(ukg^2-b^2).*tan((a/2)*sqrt(ukg^2-b^2))-sqrt(b^2+uk^2);
        bfws = vpasolve(R == 0, b, [uk0+eps, ukg-eps]);
        if isempty(bfws)
            betaFWs(i,j) = NaN;
        else
            betaFWs(i,j) = bfws(1);
        end
        %TERMINOS EQ TRASENDENTAL BETA ANTI SIMETRICO
        I(b) = sqrt(ukg^2-b^2).*cot((a/2)*sqrt(ukg^2-b^2))+sqrt(b^2+uk^2);
        %EQ Y SOLUCION TRASENDENTAL PARA BETA ANTI SIMETRICO
        bfwa = vpasolve(I == 0, b, [uk0+eps, ukg-eps]);
        if isempty(bfwa)
            betaFWa(i,j) = NaN;
        else
            betaFWa(i,j) = bfwa(1);
        end
    end
end
figure(1)

