%Numerik Praktikum

%Christpher Deitmers 1859196
%Julian Buttstädt 1851189

%Blatt 1
%Programmieraufgabe 2

%Teilaufgabe (i): n-te Näherung von e im Grenzprozess


function e = Euler1(n)
    
    for i = 1:length(n)
        e(1,i) = (1 + n(i)^(-1)).^(n(i));
    end
    
end

