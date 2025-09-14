%%%%%%%%%%%%Numerik Praktikum Blatt 5%%%%%%%%%%%%%%%%
%%% 
%%% Gruppe 4: Christopher Deitmers, 1859196
%%%           Robert Fladung, 1822623
%%%           Julian Buttstädt, 1851189


%%%        A12 - Adaptives Simpson-Verfahren
%%% 
%%% Eingabe:    Spline-Koeffizienten 
%%%             für ein Polynom 4. Grades in der Form
%%%             s(x)=ax^4+bx^3+cx^2+dx+e             
%%%             abstand abs zum nächsten fahrzeug
%%% 
%%% Ausgabe:    Länge des zurückzulegenden Wegs bis zum Zielfahrzeug in einem Zyklus,
%%%             wobei die Position des anderen Fahrzeugs in Polarkoordinaten 
%%%             gegeben ist. Also berechnet sich die Länge pro Zyklus
%%%             aus einem Intervall von 0 bis 2*pi.
%%%
%%% Wahl der Genauigkeit:
%%%             wähle die absolute Fehlerschranke abhängig vom aktuellen
%%%             Abstand (je kleiner der Abstand, desto kleiner die
%%%             Schranke), z.b. abs/100
%%%             der relativen Fehlerschranke genügen 1/100 
%%%             so sollten ausreichend genaue ergebnisse vorliegen, ohne
%%%             unnötig viel Rechenaufwand zu verursachen, wie wenn man z.B
%%%             die MaschinenGenauigkeit fordert
%%%             sei f


function [l] = sflength(a,b,c,d,e,abs)

    RTOL = 1/100; 
    ATOL = abs/100;
    

    f1 = @(t)(sqrt(1+(4*a*t^3+3*b*t^2+2*c*t+d)^2));
    
    l = adsim(f1,0,2*pi,ATOL,RTOL);

end
