%%% Numerik Praktikum
%%% Blatt 4
%%% Skript

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196


%%% Aufgabe 11.2
%%% Wird ein Grenzwert überschritten?

%%% Eingabe:    Messwerte y=[y1,...,yn] in Kelvin 
%%% Ausgabe:    Logical WARNING, die positiv ausgegeben wird, falls eine
%%%             Überschreitung von 100 Grad anzunehmen ist, 
%%%             negativ, falls nicht.

%%% Funktionsweise
%%%     1.  Berechne den Spline mit Aufgabe 1:
%%%             x:      gegeben in Milisekunden 
%%%             y:      Eingabe (verrauschte Messwerte)
%%%             yDelta: =1 für alle i nach Aufgabenstellung
%%%             S:      zwecks konservativer Abschätzung der Temperatur
%%%                     wähle kleinstmögliches S aus Konfidenzintervall
%%%                     ("natural values lie within confidence interval")
%%%     2.  Prüfe intervallweise, ob eine Überschreitung von 100 Grad, bzw
%%%         373.15 anzunehmen ist
%%%             Dabei prüfen wir pro Intervall 101 Zwischenwerte (x Werte
%%%             doppelt), also "alle 0.1 Milisekunden". Theoretisch könnte
%%%             zwischen den Stellen eine Überschreitung stattfinden, welche wir 
%%%             vernachlässigen.


function WARNING = hundreddegree(y)

    x=[0:10:1000]; %%in Milisekunden
    n=101;
    WARNING=0;
    %%% es gibt 1000/10 +1, also 101 Messwerte
    %%% und das Konfidenzintervall [101-sqrt(202),101+sqrt(202)]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S=101-sqrt(202);  % Temperatur steigt und fällt nicht linear, wir wollen Kurven
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:101
        yDelta(i)=1; %%% in Kelvin
    end
    
    [~,~,~,~,s]=smsp(x,y,yDelta,S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% jetzt wollen wir nach Grenzüberschreitungen suchen
    %%% 100 Grad Celsius sind 373,15 Kelvin
    for i=1:n-1
        x1=x(i);
        x2=x(i+1);
        xq=linspace(x1,x2);
        if any( ppval(s,xq) >= 373.15 )
            WARNING = 1;
            disp('Cool Down!')
            return
        end
    end
end

