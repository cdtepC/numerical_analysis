%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189

%%% Skript zu Aufgabe 4
%%% LR - Zerlegung


%%% Teilaufgabe (i)

A = [ 1 3 2; 4 2 9; 5 3 7];
[L,R] = LR_Pivoting0(A)


%%% Teilaufgabe (ii)

A = [1 -3 2; 5 1 8; 10 2 9];
[L,R,P,Q] = LR(A,1)

