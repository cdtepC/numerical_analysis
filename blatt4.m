%%% Numerik Praktikum
%%% Blatt 3
%%% Skript

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Teilaufgabe 11.1

%%% die Probe für Spline ohne smoothing geht auf
x=[1 7 12 15 19];
y=[8 10 7 8 7];
yDelta=[0 0 0 0 0];
S=0;
s=smsp(x,y,yDelta,S);
%%% der Quelle 
%%% http://www.tm-mathe.de/Themen/html/funnatsplines.html
%%% zufolge stimmt unser Ergebnis
%%% Achtung: auf der Seite hat sich ein Fehler eingeschliche:
%%%     in der Formulierung von s_1 wurde d(1) mit falschem Vorzeichen
%%%     versehen
%%%     in der dort vorangehenden Tabelle ist jedoch ersichtlich, dass
%%%     unser Vorzeichen das Richtige ist (sowie durch das Bild)


%%% Aufgabe 11.2

a = 340;
b = 385;
y = (b-a).*rand(1,101)+a;
WARNING = hundreddegree(y)


%%% Aufgabe 11.3

%%% import data auto.mat und wähle x=x', y=y'
%%% die erreichte ausgabe scheint sinnvoll


