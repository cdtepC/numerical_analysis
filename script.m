%%%%%%%%%%%%Numerik Praktikum Blatt 5%%%%%%%%%%%%%%%%
%%% 
%%% Gruppe 4: Christopher Deitmers, 1859196
%%%           Robert Fladung, 1822623
%%%           Julian Buttstädt, 1851189

%%% Aufgabe 12
f = @(x)(sin(x)^2);
a = 0;
b = 2*pi;
ATOL = 0.000001;
RTOL = 0.00000001;
I = adsim(f,a,b,ATOL,RTOL)
%%% analytisches Ergebnis: pi



%%% Aufgabe 13
%%% wähle zum Testen die konstante 1-Funktion, s.d. 2*pi perfekt sind
a = 0;
b = 0;
c = 0;
d = 0;
e = 1;
abs = 10;
l = sflength(a,b,c,d,e,abs)


%%% Aufgabe 14
%%% wähle yD mit \lambda = 1
%%% lineare konvergenz -> Verdopplung von n verdoppelt Genauigkeit des
%%% (explizites eulerverfahren hat konvergenzordnung 1)
%%%
yD=@(y)(y);
y0 = 1;     %%% analytisches Ergebnis: e^x, 
            %%% daher enthält L als letzten Wert e^1=e
a = 0;
b = 1;
n = 50000;

L = xeuler(a,b,n,y0,yD)





