%Numerik Praktikum

%Christpher Deitmers 1859196
%Julian Buttstädt 1851189

%Blatt 1
%Programmieraufgabe 2

%Teilaufgabe (iii)


function e = Euler2(n)

    e = 0;
    
    for i = 0:n
        e = e + (1/factorial(i));
    end
    
end
        