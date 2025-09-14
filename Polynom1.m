%Numerik Praktikum

%Christpher Deitmers 1859196
%Julian Buttstädt 1851189

%Blatt 1
%Programmieraufgabe 3

%Teilaufgabe (i): Polynomauswertung


function p = Polynom1(x, a) 

    n = length(a);
    p = 0;
    
%falls sum(a1:aN) verboten ist
    for i = 1:n
        p = p + a(i)*(x^(i-1));
    end
 
end
