%%%%%%%%%%%%Numerik Praktikum Blatt 5%%%%%%%%%%%%%%%%
%%% 
%%% Gruppe 4: Christopher Deitmers, 1859196
%%%           Robert Fladung, 1822623
%%%           Julian Buttstädt, 1851189


%%%        A12 - Adaptives Simpson-Verfahren
%%% 
%%% Eingabe:    zu integrierende Funktion f (function handle)
%%%             a,b Integral-Grenzen
%%%             ATOL absoluter Fehler
%%%             RTOL relativer Fehler 
%%%             sei I_s eine grobe Schätzung des zu berechnenden Integrals,
%%%             dann soll: 
%%%      |J-Jhut| < ((b_lokal-a_lokal)/(b_global-b_lokal))*max(15,RTOL*I_s)
%%%             dazu wähle I_s als Jhut 
%%% Ausgabe:    Approximation I des Integrals von f auf [a,b]
%%%  

function [I] = adsim(f,a,b,ATOL,RTOL)

                    %%% Berechne J(f) über die Fassregel und schätze den
                    
J = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));

                   %%% Berechnung der Jhut(f) über die Fassregel
                   %%% erfolgt in einem Schritt->Aufwand 
                   %%% Funktion soll so selten wie mögl. ausgewertet werden
Jhut = ((b-a)/12)*(f(a)+2*f((a+b)/2)+4*f((3*a+b)/4)+4*f((a+3*b)/4)+f(b));  


    if (abs(J-Jhut) < 15*ATOL) && (abs(J-Jhut) < abs(Jhut)*RTOL) 
        %%% Jhut*RTOL, um die Abweichung anteilmäßig zu beschränken
        %%% deswegen können wir RTOL im rekursiven Aufruf unverändert
        %%% übergeben
        I = Jhut;
    else
        I = adsim(f,a,(a+b)/2,ATOL/2,RTOL) + adsim(f,(a+b)/2,b,ATOL/2,RTOL);
        %%% übergebe ATOL/2 wegen Skalierung:
        %%% ((b_lokal-a_lokal)/(b_global-a_global) = (1/2)
    end
end

        