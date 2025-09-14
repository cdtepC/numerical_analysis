%Teilaufgabe 3


function y = Automat( x, f1, f2)

    for i = 1:10
        f1 = @()(f1(x^i))
        f2 = @()(f2(x^i))
        t1 = timeit(f1)
        t2 = timeit(f2)
    end
end
    
    