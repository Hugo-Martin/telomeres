function [phat, pci, h, pval, st] = test_geom(values, min_attendu, histobool, nbreal,Emin)

% Fonction qui teste si un jeu de données suit une loi géométrique.

% S'il y a trop peu de classes attendues (i.e. dans le test d'adéquation, moins de trois classes), on renvoie la valeur 0. 

% OUTPUT:
%
% phat : paramètre estimé par EMV
% pci : intervalle de confiance à 95% pour le paramètre estimé phat
% pval: valeur p renvoyée au test

% INPUT:
%
% values : liste de valeurs
% min_attendu : 0 ou 1, selon la valeur minimale théorique du support de la
% distribution, 0 si c'est une géométrique "à la matlab" ou 1 si c'est "à
% la française"
% histobool : 0 si on ne veut pas de représentation graphique, 1 si oui



if min_attendu == 1
    shift = 1;
elseif min_attendu == 0
    shift = 0;
else
    disp('Wrong minnimal expected value')
end

% values <- values - shift : dans MATLAB, une loi géométrique commence à 0

[phat, pci] = mle(values - shift,'distribution','geo'); % estimateur de maximum de vraissemblance

histcnts = histcounts(values - shift,'BinMethod','integers'); % compte du nombre de séquences de cycles longs de chaque longueur


pd = makedist('NegativeBinomial','R',1,'p',phat); % fabrication d'une distribution géométrique de paramètre estimé

bins = 0:(numel(histcnts)-1);

expected = sum(histcnts) * pdf(pd,bins);

[h,pval,st] = chi2gof(bins,'Ctrs',bins,'Frequency',histcnts,'Expected',expected,'NParams',1,'EMin',Emin);

% %% Valeurs attendues
%     
%     nb_classes = min(floor(1 + (log(5) - log(sum(histcnts)))/log(1 - phat)),numel(histcnts)); % pour un test du chi2, on dit communément qu'on doit attendre au moins 5 observations par catégories
%     expected = (1-phat).^(0:(nb_classes-1));
%     expected(1:(end-1)) = phat.*expected(1:(end-1));
%     expected = numel(values).*expected;
%     
% %% Observations
% 
%     if nb_classes > 2
%         observed = histcnts(1:nb_classes);
%         observed(end) = observed(end) + sum(histcnts((nb_classes + 1):end));
%     end
% 
%     
%     
% %% Détermination de la valeur p
%     
%     if nb_classes > 2
%         chi2 = sum(((expected-observed).^2)./expected);
%         pval = 1 - gammainc(chi2/2, (numel(expected)-2)/2); % numel(expected) classes et paramètre estimé donne numel(observed) - 2 degrés de liberté
%     else
%         pval = 0;
%     end
    
%% Représentation graphique

    if histobool == 1 % si on veut une représentation graphique
        %rand = geornd(phat,1,numel(values)); %simulations
        
        %pd_durees_seqL=makedist('NegativeBinomial','R',1,'p',phat);
       
        
        M = max(values+shift);
        
        tableau = zeros(nbreal,M+1-shift);
        exceptions = [];
        
        for k = 1:nbreal
            rand = geornd(phat,1,numel(values))+1; % +1 car en matlab une géométrique commence à 0
            r = histcounts(rand);
            if numel(r)<M+1-shift
                r = padarray(r',M+1-shift-numel(r),'post')'; % on adapte la taille de h
                tableau(k,:) = r;
            elseif numel(r)>M+1-shift
                tableau(k,:) = r(1:(M+1-shift));
                exceptions = [exceptions, rand(rand>(M+1-shift))]; % on stocke les réalisations plus grandes que le max des données expérimentales
            else
                tableau(k,:) = r;
            end
        end
        
        Q = quantile(tableau,[0.025 0.975],1);
        
        
        figure;
        histogram(values,'BinMethod','integers');
        hold on
        plot(shift:size(Q,2)-1+shift,Q(2,:),'*',shift:size(Q,2)-1+shift,Q(1,:),'*')
        %xlabel('Number of divisions','FontSize',18)
        %ylabel('Number of ')
        %legend('Experimental data','Lower bound of the 95% quantile','Upper bound of the 95% quantile')
        %title('Number of divisions between the two first long cycles of type A lineages','FontSize',24)
        hold off
        
    end