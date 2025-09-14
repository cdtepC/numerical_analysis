%Numerik Praktikum

%Christpher Deitmers 1859196
%Julian Buttstädt 1851189

%Blatt 1
%Programmieraufgabe 3

%Teilaufgabe (ii): Horner Schema Implementierung


function P = Polynom2(x, a )
    
    n = length(a)-1;
    P = a(n+1);
%a(n+1) entspricht a(n) in der Polynomdarstellung
    if n >= 1
        for i = 1:n
            P = P*x + a(n-i+1); 
        end
        
end