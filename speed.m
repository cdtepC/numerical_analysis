%%% Numerik Praktikum
%%% Blatt 4
%%% Skript

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196


%%% Aufgabe 11.3
%%% Fahrprofil eines vorbeifahrenden Fahrzeugs

%%% Eingabe: Abstands-Messdatenreihe (x,y)
%%%             y in Metern
%%%             x in Sekunden
%%% Ausgabe: Geschwindigkeitsprofil des betrachteten Objekts
%%%             in Form eines piecewise polynoms v
%%%             v in Metern pro Sekunde
%%% 

%%% Funktionsweise
%%%     Wahl eines natürliches S, wie im Skript angegeben, liefert eine
%%%     angemessene Glättung (aus dem Konfidenzintervall, konservativ)
%%%     berechne den smoothing spline zur gegebenen Messreihe und
%%%     Glättungsparameter S und festen Messungenauigkeiten von 20m
%%%     
%%%     von den nun geglätteten Informationen berechnet sich die
%%%     Geschwindigkeit als Änderungsrate der Abstände (1. Ableitung)
%%%     da es sich um ein intervallweises polynom hadelt, können wir auch
%%%     nur intervallweise ableiten

function v = speed(x,y)

    [~,n] = size(x);
    yDelta=zeros(1,n);
    for i=1:n
        yDelta(i)=20;
    end
    S= n+sqrt(2*n); %%% n-sqrt(2n) <= S <= n+sqrt(2n)
    
    abs = smsp(x,y,yDelta,S);
    plot(x,ppval(x,abs))
    hold on
    y=zeros(n-1,3);
    v = mkpp(x,y);
    
    for i=1:n-1
        for j=1:3
            v.coefs(i,j)= (j+1)*abs.coefs(i,j+1);
        end
    end
    
    plot(x,ppval(x,v));
    
end
    
    
    
    
        