%%% Numerik Praktikum
%%% Blatt 3
%%% Augabe 8

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196

%%% Eingabe: Matrix A mit A: R^n -> R^m 
%%% Ausgabe: Moore-Menrose-Pseudoinverse A_plus
%%% Funktionsweise:
%%%     führe Singulärwertzerlegung aus und erhalte A = U[Sigma;N]V'
%%%     berechne Sigma_plus und daraus A_plus




function A_plus = MPPi(A)

    [U,S,V] = svd(A);
    [m,n] = size(A);
    
    N = zeros(m-n,n);
    S_plus=zeros(n,n);
    
    %%% passe U und V an Aufgabenstellung an
    %%% U soll aus R^m->R^m, ist aber aus R^n->R^m
    %%% V soll aus R 
    
    %%% falls Eigenwerte sortiert worden sein sollen tue dies an dieser
    %%% Stelle
    
    %%% erstelle Sigma_plus
    for i=1:n
        if S(i,i)==0
            S_plus(i,i) = 0;
        else
            S_plus(i,i) = 1/S(i,i);
        end
    end
    
    A_plus = V*[S_plus,N']*U';
end
    
    
    