%%% Numerik Praktikum
%%% Blatt 3
%%% Robert Fladung
%%% Julian Buttstädt
%%% Christopher Deitmers

%%% Aufgabe 9 - Lineare Ausgleichsrechnung
%%% erste Teilaufgabe

%%% löse Normalengleichung mittels geometrischem Zugang
%%% mit Hilfe der Cholesky-Zerlegung

%%% Eingabe: funktionaler Zusammenhang bigA: R^n -> R^m mit m>n
%%%          b \in R^m             
%%% Ausgabe: Vektor xHat \in R^n, der |Ax-b| minimiert (bzgl 2-Norm)

%%% Funktionsweise: 
%%%          xHat erfüllt die Normalengleichung (A'A)x=(A'b), ein LGS,
%%%          welches wir mittels der Cholesky-Zerlegung lösen können, falls
%%%          A vollen Spaltenrang hat, denn
%%%          A'A ist symmetrisch und semi-positiv-definit
%%%          falls rg_s(A)=n, dann sogar positiv-definit
%%%          dann führen wir eine Vorwärts- und eine Rückwärtssubstitution
%%%          durch
%%%%%%%
    %%% Wenn A keinen vollen Spaltenrang hat,
    %%% dann existiert keine eindeutig bestimmte Lösung unseres Problems.
    %%% Matlab benutzt zur Bestimmung des Rangs die Singulärwertzerlegung,
    %%% welche in ihrer numerischen Bestimmung (laut Wikpedia)
    %%% einen Aufwand von ca (4/3)*n^3 + c*n^2 + b hat.
    %%% Die Choleskyzerlegung hingegen nur ca (1/6)*n^3.
    %%% Es lohnt sich also nicht, vor der Chol-Zerl. auf vollen Rang zu
    %%% prüfen.

function xHat = geomLARcz( A, b)

    %%% Vorbereitung für das LGS AHat xHat = bHat 
    bHat = A'*b;
    [L,n] = CZ_modif(A'*A);
    y = zeros(n,1);
    
    if (L == 0)
        disp('Es existiert keine eindeutig bestimmte Lösung')
        return
    end
    
    %%% Vorwärtssubstitution L y = bHat
    y(1) = bHat(1)/L(1,1);
    for i=2:n
        y(i) = (bHat(i) - L(i,1:i-1)*y(1:i-1) )/L(i,i);
    end
    
    %%% Rückwärtssubstitution L^T x = y
    xHat(n) = y(n)/L(n,n);  %%% Diagonalelemente von L und L^T sind gleich
    for i=n-1:-1:1
        y(i) = (y(i) - L(i+1:n,i)' )/L(i,i); 
    end
    xHat = xHat';
    
end
    
    
        
        
    
    
    
    
    