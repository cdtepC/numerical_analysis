%%% Numerik Praktikum
%%% Blatt 4
%%% Skript

%%% Julian Buttstädt 1851189
%%% Robert Fladung 1822623
%%% Christopher Deitmers 1859196

%%% Aufgabe 1 - smoothing spline-interpolation


%%% Eingabe: Messreihe (x,y) als arrays x und y
%%%             es muss x(i)<x(i+1) für i=1,...,(n-1)
%%%          Unsicherheiten yDelta: R^n -> R
%%%          Glättungsparameter S aus R^(+)
%%% Ausgabe: glättender, kubischer Splineinterpolant s in Form von:
%%%          struct s, der Form 'pp' (piecewise polynom)
%%%          so entsteht der Spline-Interpolant s mit
%%%                 s(x)= d(i)*x^3+c(i)*x^2+b(i)*x+a(i)
%%%                 gegeben     x(i) <= x < x(i+1) für i=1,...,n-2
%%%                 und         x(n-1) <= x <= x(n)
%%%          


%%% Funktionsweise:
%%%
%%%     theoretische Grundlage:
%%%         nach der Interpolationsforderung erfüllt der Interpolant p:
%%%         p(x(i))=y(i), damit auch s
%%%         wenn wir Messungenauigkeiten yDelta miteinbeziehen, so müssen
%%%         wir unsere Forderung anpassen zu: 

%%%     Wir übersetzen den Algorithmus aus dem gegebenen Artikel über
%%%     smoothing spline-interpolation.
%%%     Dort werden "relative weights dy(i)" gegeben, welche bei uns als
%%%     Unsicherheiten yDelta(i) vorliegen.

%%% Bemerukung:
%%%     um den Code etwas effizienter zu gestalten könnte man anstatt all
%%%     der einzelnen IndexVerschiebung einmalig (in Schleifen/..) einen
%%%     anderen Index verschieben (um nicht jedes mal erneut +1 zu
%%%     berechnen)

function s = smsp(x,y,yDelta,S)

    %%% Vorbereitung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,n]=size(x);
    n1=1;
    n2=n;
    m1=n1+1;
    m2=n2-1;
    h=x(m1)-x(n1);
    f=(y(m1)-y(n1))/h; %%% dividerte differenz
    NotFIN=1;
    p=0;
    
        
    a=zeros(1,n-1);
    b=zeros(1,n-1);
    c=zeros(1,n-1);
    d=zeros(1,n-1);
    
    
    %%% Q^T *D ist tridiagonal
    r = zeros(1,n+2); %obere
    r1 = zeros(1,n+2); %mittlere
    r2 = zeros(1,n+2); %untere
    
    %%% wir brauchen T und Q mit T*c = (Q^T)*a, T ist tridiagonal
    t = zeros(1,n+2); % mittlere
    t1 = zeros(1,n+2); % obere und untere
    
    %%% 
    u = zeros(1,n+2);
    v = zeros(1,n+2);
    %%% Im gegebenen Programmcode beginnen die ArrayIndizes 
    %%% bei 0. Weil Matlab-Arrays mit Index 1 anfangen und
    %%% um möglichst nah am Code zu bleiben bauen wir bei jedem
    %%% Ansprechen eines Arrayeintrags aus obigen vordeklarierten
    %%% Arrays r,r1,r2,t,t1,u,v einen Shift (+1) ein 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Schleife (1): stelle T und R= Q^T *D auf
    for i=m1:m2  
        g=h; 
        h=x(i+1)-x(i);
        %%% g ist das h vom vorangegangenen Durchlauf,
        %%% h_i = g_{i+1}
        e=f; 
        f=(y(i+1)-y(i))/h;
        %%% e ist f vom vorangegangenen Durchlauf, 
        %%% die Mittlere Steigung (Dividierte Differenz)
        %%% f_i = e_{i+1}
        
        a(i)=f-e; 

        t(i+1)=2*(g+h)/3;
        t1(i+1)=h/3; %%% t_i=t_{i+1,i}=t_{i,i+1}
        %%% R = Q^T *D mit r obere diag, r1 mittlere und r2 untere
        r2(i+1)=yDelta(i-1)/g;              %%% 
        r(i+1)=yDelta(i+1)/h;               %%% 
        r1(i+1)=-(yDelta(i)/g + yDelta(i)/h); %%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%% 
    
    
    %%% Schleife (2): R^T R ?: b obere diag, c mittlere, d untere
    for i=m1:m2
        b(i)=r(i+1)^2 + r1(i+1)^2 + r2(i+1)^2; %%%mittlere??
        c(i)=r(i+1)*r1(i+2) + r1(i+1)*r2(i+2); %%%innere???
        d(i)=r(i+1)*r2(i+3);                    %%% obere????
    end
    f2=-S;
        
    %%% next iteration
    while NotFIN
        
        %%%Schleife (3):sukzessive CholeskyZerlegung von((Q^T *D)*(D*Q)+pT)
        %%%                 rational version (?)
        %%%              bei gleichzeitiger Vorwärtssubstituion mit 
        %%%              R^T *(R *u) = Q^T y
        %%%              linearer rechenaufwand wegen der vorliegenden pos. def. symm.
        %%%              Band-Matrix-Form von ((DQ)^T)(DQ)
        for i=m1:m2
            r1(i)=f*r(i);
            r2(i-1)=g*r(i-1);
            r(i+1)=1/(p*b(i)+t(i+1)-f*r1(i)-g*r2(i-1)); %%% 
            u(i+1)=a(i)-r1(i)*u(i)-r2(i-1)*u(i-1); %%% Vorwärtssubstitution
            f=p*c(i) + t1(i+1) - h*r1(i);
            g=h;
            h=d(i)*p;
        end
        
        %%% Schleife (4): bestimme u nach vorangegangener Zerlegung
        %%% (Rückwärtssubstitution)
        for i=m2:-1:m1
            u(i+1)=r(i+1)*u(i+1)-r1(i+1)*u(i+2)-r2(i+1)*u(i+3);
        end
        e=0;
        h=0;
        
        %%% Schleife (5): bestimme v 
        for i=n1:m2
            g=h;
            h=(u(i+2)-u(i+1))/(x(i+1)-x(i));
            v(i)=(h-g)*yDelta(i)^2;
            e=e+v(i+1)*(h-g);
        end
        
        g=-h*yDelta(n2)^2;
        v(n2+1)=-h*yDelta(n2)^2;
        e=e-g*h;
        g=f2;
        f2=e*p^2;   %% sind vllt g und f2 die neg. bzw pos. wurzel von S?
        
        %%% if then go to fin
        if (f2>=S)||(f2<=g)
            NotFIN=0;
            continue;
        end
        
        f=0;
        h=(v(m1+1)-v(n1+1))/(x(m1)-x(n1));
        
        %%% Schleife (6)
        for i=m1:m2
            g=h;
            h=(v(i+2)-v(i+1))/(x(i+1)-x(i)); %% dividierte differenz von v
            g=h-g-r1(i)*r(i)-r2(i-1)*r(i-1);
            f=f+r(i+1)*g^2;
            r(i+1)=g;
        end
        
        h=e-p*f;      
              
        %%% if then go to fin 
        if h<=0
            NotFIN=0;
            continue;
        end
        p=p + (S-f2)/((sqrt(S/e) + p)*h); %%% etwas anders wegen NewtonMethod
        %%% anders als der gegebene algorithmus im abschnitt
        %%% "4. determination of the lagrangian parameter"
        %%% darauf wird im skript nicht genau eingegangen und die NewtonMethod
        %%% wurde bisher nicht behandelt. vorher wird sie wohl schon
        %%% benutzt, um die positive Wuzel von S zu bestimmen 
        %%% im Skript heißt es früher man verfährt dabei mit F(1/p) ??? 
        %%%
        %%% man bestimmt p nach der vorliegenden gleichung: (Falls
        %%% F(0)>sqrt(S) )
        %%% F(p)= || DQ(Q^TD^2+pT)^{-1}Q^T y||_2 = sqrt(S)
        %%% bzw apply F(1/pHat) to the NewtonMethod
        %%% wäre dann also || DQ(Q^TD^2+(1/pHat)T)^{-1}Q^T y||_2 = sqrt(S)
        %%%             zu lösen?
    end
        
    %%% fin
    %%% mit beiden Datenpaketen und Parameter S=100, bzw 70, gibt es nur 
    %%% 2 bzw. 1 Durchlauf der while-Schleife -> sehr klein
    
    %%% Schleife (7): berechne a und c
    for i=n1:n2
        a(i)=y(i)-p*v(i+1); % aus (12)
        c(i)=u(i+1); %% aus (13)
    end
    
    %%% Schleife (8): berechne b,d nach (8),(9)
    for i=n1:m2
        h=x(i+1)-x(i);
        d(i)=(c(i+1)-c(i))/(3*h);
        b(i)=(a(i+1)-a(i))/h - (h*d(i)+c(i))*h;
    end
    
    
    %%% Erweiterung, um function handle zurückzugeben und für Ansatz zu A3
    coefs = zeros(m2,4);   
    for i=n1:m2
        coefs(i,1)=d(i);
        coefs(i,2)=c(i);
        coefs(i,3)=b(i);
        coefs(i,4)=a(i);
    end
    
    s=mkpp(x,coefs);
    xq=linspace(x(n1),x(n2));
    plot(xq,ppval(s,xq))
    
    
    
end
        
        
    
    

    

