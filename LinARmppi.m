%%% Numerik Praktikum
%%% Blatt 3
%%% Augabe 8

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196

%%% Eingabe: Matrix A mit A: R^n -> R^m 
%%%          Vektor b
%%% Ausgabe: x, der Lineares Ausgleichsproblem min|Ax-b| löst
%%% Funktionsweise:
%%%     führe Singulärwertzerlegung aus und erhalte A = U[Sigma;N]V'
%%%     berechne Sigma_plus und daraus A_plus
%%%     dann lässt sich das LAR lösen duch x = A_plus b



function x = LinARmppi(A,b)

    A_plus = MPPi(A);
    x = A_plus*b;
    
end
    
