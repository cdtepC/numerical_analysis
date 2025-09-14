%Numerik Praktikum
%Julian Buttstädt 1851189
%Christopher Deitmers 1859196
%Robert Fladung 1822623
%Aufgabe 5

%Cholesky-Zerlegung
function [L]=CZ(a)
g=size(a);
n=g(1);

for i=1:n
    l(i,i)=a(i,i);
    for k=1:i-1
        l(i,i)=l(i,i)-abs(l(i,k))^2;
    end
    l(i,i)=sqrt(l(i,i));
    
    for j=(i+1):n
        l(j,i)=a(j,i);
        for k=1:i-1
            l(j,i)=l(j,i)-(l(j,k)*l(i,k));
        end
        l(j,i)=l(j,i)/l(i,i);
    end
end

L=l;

end