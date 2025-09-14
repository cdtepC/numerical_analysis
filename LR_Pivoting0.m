%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189

%%% Aufgabe 4
%%% LR - Zerlegung
%%% Teilaufgabe (i)
%%% LR-Zerlegung ohne Pivotisierung

%%% Eingabe: Matrix A mit A: R^n -> R^n
%%% Ausgabe: linksuntere Dreiecksmatrix L
%%%          rechtsobere Dreiecksmatrix R
%%%          mit A = L*R

%%% Funktionsweise: 
%%% Anwendung der Gauß-Elimination mit Speicherung der für die 
%%% Eliminationsschritte verwendeten, invertierten Elementarmatrizen, um L
%%% zu bestimmen
%%% R ergibt sich direkt aus der Anwendung auf A

function [L,R] = LR_Pivoting0(A)

    [m,n] = size(A);
    L = eye(n);
    R = A;
    
    for j = 1:n-1         
        for k = j+1:n              
            %%% bestimme Elementarmatrix N(-(R_kj/R_jj))
            N = eye(n);
            N(k,j) = -( R(k,j)/R(j,j)); %%%R(k,j)/R(i,j) = l_ij 
            %%% passe L an
            L(k,j) = R(k,j)/R(j,j);
            %%% eliminiere Eintrag R(j,k)
            R = N*R; 
        end
    end
    
end

    
