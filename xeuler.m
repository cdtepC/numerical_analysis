%%%%%%%%%%%%Numerik Praktikum Blatt 5%%%%%%%%%%%%%%%%
%%% 
%%% Gruppe 4: Christopher Deitmers, 1859196
%%%           Robert Fladung, 1822623
%%%           Julian Buttstädt, 1851189


%%%        A14 - explizites Euler-Verfahren
%%% 
%%% Eingabe:    Intervall [a,b] bzw. die einzelnen Grenzen a,b
%%%             n Anzahl der Diskretisierungspunkte 
%%%             DGL-Parameter:      yD als function-handle
%%%                                 y_0 als Anfangswert
%%% 
%%% Ausgabe:    Lösungsapproximation L als 2xn Vektor mit 
%%%                         Zeile 1: Diskretisierungspkt. t_i
%%%                         Zeile 2: y(t_i) (y numerische Lösung der DGL)
%%%                                  ausgewertet auf den Diskr.-Pkt. t_i

function [L]=xeuler(a,b,n,y0,yD) 

    h=(b-a)/n; %%% Schrittweite
    L = zeros (2,n);
    
    %%% Testbedingungen (DGL):y'=y führt analytisch auf die Exponentialfkt.
    %%%                                         (unter Anwendung von y0=1)

    
    L(1,1) = a;
    L(2,1) = y0;
    for i = 1:n-1
        L(1,i+1) = a+i*h; % Diskretisierungsknoten (Stellen)
        L(2,i+1) = L(2,i)+h*yD(L(2,i)); % (approx. Funktionswerte)
        
    end
    
end