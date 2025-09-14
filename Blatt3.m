%%% Numerik Praktikum
%%% Blatt 3
%%% Skript

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aufgabe 7
A = [1 1 1; 0.01 0 0.01; 0 0.01 0.01];
[Q,R] = GR(A)
%%% läuft für quadratische
A = [3 5;0 2;0 0;4 5];
[Q,R] = GR(A);
%%% und auch für nicht quadratische


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aufgabe 8
%%% berechne Moore-Menrose-Pseudoinverse A_plus
A = [1 2;3 4;5 6;7 8];
A_plus = MPPi(A);
%%% läuft
A = [1 0 0;-4 1 0;0 0 1];
A_plus = MPPi(A);
%%% läuft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aufgabe 9
%%% Löse Lineare Ausgleichsrechnung 
%%%
A = [1 2;3 4;5 6;7 8];
b = [1;3;3;5];
%%% soll: [0;3/5] 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Teilaufgabe 1
%%% via A'Ax=A'b mit Cholesky
x = LinARcz(A,b)
%%% passt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Teilaufgabe 2
%%% via Givens Rotation
x = LinARgr(A,b)
%%% passt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Teilaufgabe 3
%%% via x = A_plus b mit Moore-Penrose-Pseudoinversen aus Aufgabe 8
x = LinARmppi(A,b)
%%% passt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aufgabe 10

