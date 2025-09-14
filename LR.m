%%% Numerik Praktikum
%%% Blatt 2
%%% Christopher Deitmers 1859196
%%% Julian Buttstädt 1851189

%%% Aufgabe 4
%%% LR - Zerlegung
%%% Teilaufgabe (ii)
%%% LR-Zerlegung mit möglicher Pivotisierung

%%% Eingabe: Matrix A mit A: R^n -> R^n
%%%          varargin als LogicalVariable pivoting
%%% Ausgabe: linksuntere Dreiecksmatrix L
%%%          rechtsobere Dreiecksmatrix R
%%%     falls pivoting = 0
%%%          mit A = L*R
%%%     falls pivoting = 1 
%%%          P, die Zeilenpermutaionsmatrix
%%%          mit P*A = L*R
%%%     falls pivoting = 2
%%%          P, die Zeilenpermutaionsmatrix
%%%          Q, die Spaltenpermutationsmatrix
%%%          mit P*A*Q = L*R

%%% Funktionsweise:
%%% Falls pivoting = 0, benutze LR_Pivoting0, sowie wenn keine optionale
%%% Angabe gemacht wird.
%%% Falls pivoting = 1 benutze LR_PivotingTotal.
%%% Falls pivoting = 2 benutze LR_PivotingColumn.


function [L,R,P,Q] = LR(A,varargin) 
%%% varargin bietet optionale Eingabeparameter, hier für die Wahl des Pivoting

    pivoting = varargin{1};

    %%% falls keine oder falsche Angabe zum Pivoting gemacht wird
    if nargin < 2
        pivoting = 0;
    elseif pivoting > 2 || pivoting < 0
        disp("falsche Eingabe")
        return;
    end    
    
    switch pivoting   
        %%% ohne pivoting 
        case 0
            [L,R] = LR_Pivoting0(A);
            P = [];    
            Q = [];
        %%% mit totaler pivotisierung    
        case 1 
            [L,R,P,Q] = LR_PivotingTotal(A);
        %%% mit Spaltenpivotisierung
        case 2
            [L,R,P] = LR_PivotinColumns(A);
            Q = [];
    end
end


    
    
    
    
    
    