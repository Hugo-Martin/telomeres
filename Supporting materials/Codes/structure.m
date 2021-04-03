function [seq,type] = structure(div, criterion)

%fonction qui determine la "structure" d'une lignée : nombre de cycles dans
%chacun des états (longs ou courts) successifs, séquences ordonnées
%OUTPUT
%seq : un vecteur composé des longueurs en nombre de générations des
%séquences de cycles longs ou courts
%type : vecteur de même taille composé de 0 (séquence de courts) et de 1
%(séquence de longs)
%INPUT
%div : vecteur des longueurs des cycles d'une lignée
%criterion : critère pour départager cycle long/court

bool = div>criterion;
state = bool(1); %état actuel (cycle long ou court)

seq = [];
type = state;

old = 1; %compteur pour les séquences
new = 1;



while isempty(new) == 0 %tant qu'on n'a pas atteint la fin du cycle
    new = find(bool(old:end) ~= state,1); %on trouve le premier indice d'une état différent
    if isempty(new) == 0
        seq = [seq; new-1]; %remplissage de seq
        old = old + new -1; %mise à jour de l'indice
        state = bool(old); %mise à jour de l'état
        type = [type; state]; %remplissage de type
    end
end
%disp(old)
seq = [seq; numel(div)-old+1]; %dernier remplissage de seq


