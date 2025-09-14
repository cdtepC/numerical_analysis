%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623

%%% Aufgabe 4
%%% LR - Zerlegung
%%% Teilaufgabe (ii) - Unterfunktion
%%% LR-Zerlegung mit totaler Pivotisierung: P*A*Q = L*R

%%% Eingabe: Matrix A mit A: R^n -> R^n
%%% Ausgabe: linksuntere Dreiecksmatrix L
%%%          rechtsobere Dreiecksmatrix R
%%%          P, die Zeilenpermutaionsmatrix
%%%          Q, die Spaltenpermutationsmatrix

%%% Funktionsweise:
%%% Speichere in jedem Schritt k 
%%% das Produkt N_k  aus n-k Elementarmatrizen derselben Art
%%% und die Permutationsmatrix P_k, die sich aus 2 Transpositionsmatrizen
%%% zusammensetzt. 
%%% Danach berechne L aus P_1 bis P_k und N_1 bis N_(n-k),
%%% R hat sich dann bereits analog zu dem R aus LR_Pivot0 ergeben. 

%%% Hinweis zur Notation:
%%% Im Skript beschreibt N eine Matrix, die eine Zeilenaddition durchführt. 
%%% Im Skript im Beweis von Satz 3.2 werden als K, Produkte solcher N
%%% bezeichnet.
%%% Da in unserem Programm die Speicherung der einzelnen
%%% Zeilenadditionsmatrizen nicht wichtig ist, werden bei uns die Produkte
%%% als N bezeichnet (Skript K -> N Programm)
%%% Desweiteren werden dort im Skript mit K\tilde Produkte 
%%% aus K(Skript) und Transpositionsmatrizen bezeichnet.
%%% In unserem Programm bezeichnen K die invertierten K\tilde aus dem Skript
%%% (Skript K\tilde -> K^(-1) Programm).


function [L,R,p,Q] = LR_PivotingTotal(A)

    [~,n] = size(A);
    L = eye(n);
    R = A;
    
    %%% wir brauchen n-1 Schritte, also n-1 Produkte von Elementarmatrizen
    %%% und n-1 Permutationsmatrizen, also (n-1)*2 Transpositionsmatrizen

    P = zeros(n,n, n-1, 2);
    N = zeros(n,n, n-1); %%% jede Matrix N_k ist wiederum Produkt aus n-k
                        %%% Elementarmatrizen derselben Art, welche nicht
                        %%% extra gespeichert werden

                        
                        
    %%% Zerlegung

    RR = R;
    %%% bestimme R, alle N und alle P
    for k = 1:n-1         
               
        %%% Pivotisierung 
        
         
        
        %%% bestimme max Eintrag von R
                %%% hole max aus restlicher Teilmatrix, ohne Zeilen/Spalten
                %%% 1,...,k-1. Das verursacht eine Verschiebung der Indizes
            if k > 1
                RR = R;
                for y = 1:k-1
                    RR(y,:) = [];
                    RR(:,y) = [];
                end
            end
            
        [maxRz, maxRIndexV] = max(RR);
        %%% maxRz liefert spaltenweise maxima
        [pivot, pivotIndex2] = max(maxRz);
        pivotIndex1 = maxRIndexV(pivotIndex2)+k-1;
        pivotIndex2 = pivotIndex2 + k-1;

        if pivot == 0
            disp("Eingabe unzulässig")
            return;
        end
        
            
        
        
        %%% bestimme P_k 
        %%% P_k1 Zeilentransposition 
        P(:,:,k,1) = eye(n); 
        P(pivotIndex1,pivotIndex1,k,1) = 0;
        P(k,k,k,1) = 0;
        P(pivotIndex1,k,k,1) = 1;
        P(k,pivotIndex1,k,1) = 1;       
        %%% P_k2 Spaltentransposition
        P(:,:,k,2) = eye(n);
        P(k,k,k,2) = 0;
        P(pivotIndex2,pivotIndex2,k,2) = 0;
        P(pivotIndex2,k,k,2) = 1;
        P(k,pivotIndex2,k,2) = 1;       
        
        %%% permutiere A, sodass max(A) zum Pivotelement wird
        R = P(:,:,k,1)*R*P(:,:,k,2);
        
        
        
        %%% Elimination der Einträge unter dem Pivotelement 
        
        %%% bestimme N_k
        N(:,:,k) = eye(n);
        for i = k+1:n
            N(i,k,k) = -( R(i,k)/pivot );          
        end

            %%% Hilfestellung %%%%%%%%%%%%%%%%%%%%
            if pivot ~= R(k,k)
                disp("Programm kaputt")
                return;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        
        %%% Elimination der Einträge unter (k,k)
        R = N(:,:,k)*R;
        
    end
    
    
    
    
    
    %%% bestimme L aus K(1) bis K(n-1), die aus P und N entstehen
    
    %%% bestimme K 
    K = zeros(n,n, n-1);
    K(:,:,n-1) = inv(N(:,:,n-1));  %%%invertiere N
    for i = 1:n-2
        K(:,:,i) = (-1).*N(:,:,i) + 2.*eye(n);
        for z = i+1:n-1
            K(:,:,i) = P(:,:,z,1)*K(:,:,i)*P(:,:,z,1);   %%% P = P^(-1)
        end
    end
               
    %%% setze L aus K(1) bis K(n-1) zusammen
    L = K(:,:,n-1);
    for i = 2:n-1
        L = K(:,:,n-i)*L;
    end
        
    
    %%% bestimme P und Q
    p = eye(n);
    Q = eye(n);
    for k = 1:n-1
        p = p*P(:,:,k,1);
        Q = Q*P(:,:,k,2);
    end
    
 
     
    
end
    
    
    
    
    