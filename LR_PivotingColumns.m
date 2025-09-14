%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623

%%% Aufgabe 4
%%% LR - Zerlegung
%%% Teilaufgabe (ii) - Unterfunktion
%%% LR-Zerlegung mit Spaltenpivotisierung: P*A = L*R

%%% Eingabe: Matrix A mit A: R^n -> R^n
%%% Ausgabe: linksuntere Dreiecksmatrix L
%%%          rechtsobere Dreiecksmatrix R
%%%          P, die Zeilenpermutaionsmatrix

%%% Funktionsweise:
%%% Speichere in jedem Schritt k 
%%% das Produkt N_k  aus n-k Elementarmatrizen derselben Art
%%% und die Spaltentranspositionsmatrix P_k 
%%% Danach berechne L aus P_1 bis P_k und N_1 bis N_(n-k),
%%% R hat sich dann bereits analog zu dem R aus LR_Pivot0 ergeben. 

%%% Hinweis zur Notation: analog zu LR_PivotingTotal



function [L,R,p] = LR_PivotingColumns(A)

    [~,n] = size(A);
    L = eye(n);
    R = A;
    
    %%% wir brauchen n-1 Schritte, also n-1 Produkte von Elementarmatrizen
    %%% und n-1 Permutationsmatrizen, also (n-1)*2 Transpositionsmatrizen

    P = zeros(n,n, n-1);
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
        pivot = maxRz(1);

        pivotIndex1 = maxRIndexV(1)+k-1; %%% +k-1 weil schon k-1 Zeilen gestrichen wurden
        %%% pivotIndex1 Zeilenindex des Pivotelements im k-ten Schritt
        if pivot == 0
            disp("Eingabe unzulässig")
            return;
        end
        
        
        %%% P_k Zeilentransposition
        P(:,:,k) = eye(n);
        P(k,k,k) = 0;
        P(pivotIndex1, pivotIndex1,k) = 0;
        P(pivotIndex1,k,k) = 1;
        P(k,pivotIndex1,k) = 1;
        
        
        %%% permutiere A, sodass max(Spalte von A' bzw R) zum Pivotelement wird
        R = P(:,:,k)*R;
        
        
        
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
    
    
    
    
    
    %%% bestimme L aus K(1) bis K(n-1), die aus P_k1 und N entstehen
    
    %%% bestimme K 
    K = zeros(n,n, n-1);
    K(:,:,n-1) = inv(N(:,:,n-1));  %%%invertiere N
    for i = 1:n-2
        K(:,:,i) = (-1).*N(:,:,i) + 2.*eye(n);
        for z = i+1:n-1
            K(:,:,i) = P(:,:,z)*K(:,:,i)*P(:,:,z);   %%% P = P^(-1)
        end
    end
               
    %%% setze L aus K(1) bis K(n-1) zusammen
    L = K(:,:,n-1);
    for i = 2:n-1
        L = K(:,:,n-i)*L;
    end
        
    
    %%% bestimme P 
    p = eye(n);
    for k = 1:n-1
        p = p*P(:,:,k);
    end
    
 
     
    
end
    
    
    
    
    