%Calcule l'erreur en norme 2 moyenne commise pour un couple de paramètres
%(a,L_min) donné.

function e=erreur_moyenne(L_min,repet,a,sorted_senescence_exp,support,densite,repartition,modele,modif_sigma)


% INPUT:
%L_min : seuil de la sénescence
% repet : nombre sur lequel on va faire la moyenne
% a : paramètre dans L1 + a*L2
% beta : paramètre du modèle
% L0 : idem
% sorted_senescence_exp : liste triée des valeurs expérimentales
% d'apparition de sénescence dans les lignées, vecteur ligne

% OUTPUT:
% e : erreur associée au couple de paramètres (a,L_min)

%tic;

%de=[10 11 15 16 20 24 25 26 26 27 29 30 31 32 34 35 36 38 41 41 41 42 42 43]; % donnees expérimentales
%nb=length(de); %nombre de cellules étudiées par expérience

nn=zeros(numel(sorted_senescence_exp)*repet,1);

for k=1:numel(sorted_senescence_exp)*repet

nn(k)=generation_senescence(L_min,a,support,densite,repartition,modele,modif_sigma);

end


%nn=sort(nn,1);
nn=sort(reshape(nn,repet,numel(sorted_senescence_exp)),2);
repetde = repmat(sorted_senescence_exp,repet,1);

e = norm(reshape(nn-repetde,numel(nn),1),2)/(sqrt(repet)*norm(sorted_senescence_exp,2));
%e=1/repet*norm(nn-repmat(sorted_senescence_exp',1,repet),2);

%e=norm(reshape(nn-repmat(sorted_senescence_exp',1,repet),numel(nn),1),2)/(sqrt(repet)*norm(sorted_senescence_exp,2));



%e = norm(sort(nn,2) - reshape(repmat(sorted_senescence_exp,repet,1),repet*numel(sorted_senescence_exp),1),2)/(sqrt(repet)*norm(sorted_senescence_exp,2));

%toc;
