% Goal: studying the robustness of the description of the number of
% consecutive arrests by a geometric law, either fort internal sequences or
% terminal sequences


% a lineage is said senescent if it ends by a death (ended(j) = 1) and its
% last cycle if long (v(end) > threshold)

close all

fig_properties

min_length_traj=3;

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;

addpath('./Supporting material/Codes')
addpath('./Supporting material/Data')

%ended_exp=0; %1 si on ne veut que les lignées terminées, 0 sinon
DOX_exp=1; %1 si on veut après l'addition de DOX, 0 si on prend tout.

[cycle_lengths, indeces, ~, ended]=extraction(data_exp,DOX_exp,min_length_traj);


%
nbreal = 100;

tested_threshold = 10:35;

display_individual_plots = 1;

vecteur_p = zeros(1,numel(tested_threshold));
vecteur_pEMV = zeros(1,numel(tested_threshold));

vecteur_p_sene = zeros(1,numel(tested_threshold));
vecteur_pEMV_sene = zeros(1,numel(tested_threshold));

%% Préparation des données

for i = 1:numel(tested_threshold)
    
    threshold = tested_threshold(i);


    durees_seqL = [];
    durees_senescence = [];

    for j = 1:numel(data_exp)
        v = cycle_lengths(indeces==j);
        [s,t] = structure(v,threshold);
        t_all=t;
        if t(end) == 1 && ended(j) == 1 % si la lignée est sénescente
            durees_senescence = [durees_senescence, s(end)]; % On stocke les durées des sénescences (en générations)
            s = s(1:(end-1)); % on retire la sénescence
            t = t(1:(end-1)); % idem
            
        end
          durees_seqL = [durees_seqL; s(t == 1)]; % on stocke les longueurs des séquences des cycles longs
    end
    %disp(seuil)
    
[vecteur_pEMV(i),~ , ~, vecteur_p(i), ~] = test_geom(durees_seqL, 1, 0, 100, 5);

[vecteur_pEMV_sene(i), ~, ~, vecteur_p_sene(i), ~] = test_geom(durees_senescence, 1, 0, 100, 5);

    
    table_counts_internal = [];
    table_counts_sene = [];
    
    for k = 1:nbreal
        
        rand = geornd(vecteur_pEMV(i),1,numel(durees_seqL))+1; % +1 since geometric laws in matlab start at 0
 
        h = histcounts(rand,'BinMethod','integers');

        if numel(h) < size(table_counts_internal,2)
            h = [h, zeros(1,size(table_counts_internal,2)-numel(h))];
            table_counts_internal = [table_counts_internal;h];
        elseif numel(h) > size(table_counts_internal,2)
            table_counts_internal = [table_counts_internal, zeros(size(table_counts_internal,1),numel(h) - size(table_counts_internal,2))];
            table_counts_internal = [table_counts_internal;h];
        else
            table_counts_internal = [table_counts_internal;h];
        end
        
        
        rand_sene = geornd(vecteur_pEMV_sene(i),1,numel(durees_senescence))+1;
        h_sene = histcounts(rand_sene);
        if numel(h_sene) < size(table_counts_sene,2)
            h_sene = [h_sene, zeros(1,size(table_counts_sene,2)-numel(h_sene))];
            table_counts_sene = [table_counts_sene;h_sene];
        elseif numel(h_sene) > size(table_counts_sene,2)
            table_counts_sene = [table_counts_sene, zeros(size(table_counts_sene,1),numel(h_sene) - size(table_counts_sene,2))];
            table_counts_sene = [table_counts_sene;h_sene];
        else
            table_counts_sene = [table_counts_sene;h_sene];
        end
    end
    
    Q_durees_seqL = quantile(table_counts_internal,[0.025 0.975],1);
    Q_durees_sene = quantile(table_counts_sene,[0.025 0.975],1);
    
 
if threshold == 18
    figure;
    subplot(2,1,1)
    histo = histogram(durees_seqL,'BinMethod','integers');
    hold on
    plot(1:size(Q_durees_seqL,2),Q_durees_seqL(1,:),'*',1:size(Q_durees_seqL,2),Q_durees_seqL(2,:),'*',1:max(durees_seqL),numel(durees_seqL).*geopdf(0:(max(durees_seqL)-1),vecteur_pEMV(i)),'*')
    
    xlabel('Number of long cycles')
    ylabel('Number of occurences')
    legend('Observed frequency','2.5%','97.5%','Expected values')
    title('Number of consecutive internal long cycles ')
    text(0.6,0.7,['threshold = ' num2str(threshold)],'Units','normalized')
    text(0.6,0.6,['pEMV = ' num2str(vecteur_pEMV(i))],'Units','normalized')
    text(0.6,0.5,['p-value = ' num2str(vecteur_p(i))],'Units','normalized')
    hold off
    
    subplot(2,1,2)
    histo_sene = histogram(durees_senescence,'BinMethod','integers');
    hold on
    plot(1:size(Q_durees_sene,2),Q_durees_sene(1,:),'*',1:size(Q_durees_sene,2),Q_durees_sene(2,:),'*',1:max(durees_senescence),numel(durees_senescence).*geopdf(0:(max(durees_senescence)-1),vecteur_pEMV_sene(i)),'*')

    xlabel('Number of long cycles')
    ylabel('Number of occurences')
    legend('Observed frequency','2.5%','97.5%','Expected values')
    title('Number of consecutive terminal long cycles ')
    text(0.6,0.7,['seuil = ' num2str(threshold)],'Units','normalized')
    text(0.6,0.6,['pEMV = ' num2str(vecteur_pEMV_sene(i))],'Units','normalized')
    text(0.6,0.5,['p-value = ' num2str(vecteur_p_sene(i))],'Units','normalized')
    hold off
    
end

end


figure(13);
%subplot(2,1,1)
plot(tested_threshold*10,vecteur_p,'o',tested_threshold*10,0.05*ones(size(vecteur_p)),'MarkerSize',markersize,'LineWidth',linewidth,'MarkerEdgeColor',ColorFit1,'MarkerFaceColor',ColorFit1)
xlabel('Value of the threshold (min)','Fontsize',fontsize)
%ylabel('p-valeur')
%legend({'p-valeur','0.05'},'Location','northwest','Fontsize',fontsize)
%title('Internal sequences on long cycles','Fontsize',fontsize)
%print -depsc robustesse_geom_deuxieme_regime
ax_properties
xlim([min(tested_threshold*10) max(tested_threshold*10)])
ylim([0 1])
figure(21)
%subplot(2,1,2)
plot(tested_threshold*10,vecteur_p_sene,'o',tested_threshold*10,0.05*ones(size(vecteur_p_sene)),'MarkerSize',markersize,'LineWidth',linewidth,'MarkerEdgeColor',ColorFit1,'MarkerFaceColor',ColorFit1)
xlabel('Value of the threshold (min)','Fontsize',fontsize)
%ylabel('p-valeur')
%legend({'p-value','0.05'},'Location','northwest','Fontsize',fontsize)
%title('Senescence','Fontsize',16)
%print -depsc robustesse_geom_deuxieme_regime
ax_properties
savefig('fig21.fig')
figure(14);
%subplot(2,1,1)
plot(tested_threshold*10,vecteur_p_all,'o',tested_threshold*10,0.05*ones(size(vecteur_p_all)),'MarkerSize',markersize,'LineWidth',linewidth,'MarkerEdgeColor',ColorFit1,'MarkerFaceColor',ColorFit1)
xlabel('Value of the threshold (min)','Fontsize',fontsize)
%ylabel('p-valeur')
%legend({'p-valeur','0.05'},'Location','northwest','Fontsize',fontsize)
title('sequences of long cycles, internal or senescent','Fontsize',fontsize)
%print -depsc robustesse_geom_deuxieme_regime
ax_properties
xlim([min(tested_threshold*10) max(tested_threshold*10)])
ylim([0 1])


