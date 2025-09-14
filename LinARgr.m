%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623

%%% Hausaufgabe 9, Teilaufgabe 2


%%% Eingabe: Matrix A aus R^(mxn)
%%%          Vektor b aus R^m
%%% Ausgabe: Vektor xHat aus R^n,
%%%          der das gegebene LAR löst
%%% Funktionsweise:
%%%      finde Matrizen Q, R (Q orth., R obere Dr.), sd. A=Q*R, A also mittels
%%%      Givens-Rotation zerlegt ist und verändere sukzessiv b, sodass das
%%%      LAR mit einer Rückwärtssubstitution gelöst wird
%%%      Grundlage: min|Ax-b| = min|QAx-Qb| = min|Rx-Qb|
%%%                 dann x = HatR^(-1)*(Qb) 
%%%                 für R = [HatR; [0]]

%%% Programm GR kopiert mit folgenden Änderungen:
%%% Programmkopf
%%% Zeile 57: Givens-Rotation für b
%%% Zeile 65: Rückwärtssubstitution mit x = R^(-1)b 


function [x] = GR(A,b)

    [m,n] = size(A);
    l = min(m-1,n);
    R = A;
    Q = eye(m);
    for j=1:l         %Spalten
        for i=j+1:m     %Zeilen
            G = eye(m);
            
            if R(i,j) == 0  %Spezialfall: Pivot ist bereits 0, hier ist nichts zu tun.
                
            else
                %Hilfsvariablen t,v zur Bestimmung von c und s.               
                t = abs(R(j,j))+abs(R(i,j));             
                v = t*sqrt((R(j,j)/t)^2+(R(i,j)/t)^2);
                
                %Bestimmung von c, s, um G_ij zu bestimmen
                c = R(j,j)/v;      
                s = R(i,j)/v;
                
                G(j,j)=c;
                G(i,i)=c;
                G(i,j)=-s;
                G(j,i)=s;
                
                %Durchführung der Givens-Rotation
                Q = G*Q;            
                R = G*R;   
                b = G*b;                
                            
            end 
        end
       
        
    end
    %%% löse das LAR via x=HatR^(-1)*b
    Q=Q';
    x = [inv(R(1:n,:)),zeros(n,m-n)]*b;
    
end
                