% Boxplots to study how long the arrests are


close all

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

fig_properties

min_length_traj=3;

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;


%ended_exp=0; %1 si on ne veut que les lignées terminées, 0 sinon
DOX_exp=1; %1 si on veut après l'addition de DOX, 0 si on prend tout.



threshold = 18;

[cycle_lengths, indeces, ~, ended]=extraction(data_exp,DOX_exp,min_length_traj);

last_of_sene = [];
last_of_others = [];

others_of_sene = [];
others_of_others = [];

durees_exclues = [];
nb_gene_exclues = [];

for j = 1:numel(data_exp)
    v = cycle_lengths(indeces == j);
    %% Traitement à part si la lignée se termine par une séquence de longs
    if v(end) > threshold % si le dernier cycle est long
        if ended(j) == 1 % si la lignée est terminée à la fin de l'expérience, et donc sénescente
            last_of_sene = [last_of_sene, v(end)];
            v(end) = [];
            if v(end) > threshold % s'il y a encore un cycle long à la fin
                vect = flip(v);
                vect = vect(1:find(vect <= threshold, 1)-1);
                others_of_sene = [others_of_sene, vect'];
                v = v(1:end-numel(vect)); % on retire le reste de la sénescence
            end
        else % si la lignée n'est pas terminée à la fin, on ne sait pas si la séquence courante est la sénescence ou pas, donc on l'exclut
            durees_exclues = [durees_exclues, v(numel(v) - find(flip(v) <= threshold, 1) + 2:end)']; % on stocke les séquences exclues
            nb_gene_exclues = [nb_gene_exclues, find(flip(v) <= threshold, 1) - 1]; % nombre de générations dans les séquences exclues
            v(numel(v) - find(flip(v) <= threshold, 1) + 2:end) = [];
        end
    end
    %% Traitement en initialisant avec une séquence de courts à la fin
    while isempty(v) == 0 % tant qu'on n'a pas examiné toute la lignée
        i = find(v > threshold,1);
        k = find(v(i:end) <= threshold,1);
        if  isempty(k) % si on ne trouve pas de cycle court après le i-ème cycle (qui est long)
            if ~isempty(i) % si on trouve encore des cycles longs
                vect = v(i:end)';
                last_of_others = [last_of_others, vect(end)];
                others_of_others = [others_of_others, vect(1:end-1)];
            end
            v = [];
        else
            vect = v(i:i+k-2)';
            last_of_others = [last_of_others, vect(end)];
            others_of_others = [others_of_others, vect(1:end-1)];
            v = v(i+k-1:end);
        end
    end
end

figure;

%subplot(2,2,1)
boxplot([last_of_others, last_of_sene],[repmat({'Internal sequences'},1,numel(last_of_others)),repmat({'Senescence'},1,numel(last_of_sene))],'LabelVerbosity','minor')
title('Last cycle in the sequences of arrests','Fontsize',16)

figure(1920)%)%,'BoxStyle','filled',
h1=boxplot([others_of_sene, last_of_sene, others_of_others,last_of_others],[repmat({'Previous cycles'},1,numel(others_of_sene)),repmat({'Last cycles'},1,numel(last_of_sene)),repmat({'Previous cycles'},1,numel(others_of_others)),repmat({'Last cycles'},1,numel(last_of_others))],'LabelVerbosity','minor');
set(h1(:,:),'LineWidth',linewidth)
set(h1(:,:),'Color',ColorData2)
ax_properties
ylim([0 300])
savefig('Fig1920.fig')

figure(19)
h=boxplot([others_of_others, last_of_others],[repmat({'Previous cycles'},1,numel(others_of_others)),repmat({'Last cycles'},1,numel(last_of_others))],'LabelVerbosity','minor');
set(h(:,:),'LineWidth',linewidth)
set(h(:,:),'Color',ColorData2)
%title('Internal sequences','Fontsize',16)
ax_properties
ylim([0 300])
title('Internal sequences of arrests')
savefig('Fig19.fig')

figure
%subplot(2,2,3)
boxplot([others_of_others, others_of_sene],[repmat({'Other sequences'},1,numel(others_of_others)),repmat({'Senescence'},1,numel(others_of_sene))],'LabelVerbosity','minor')
title('Cycles before the last of sequences of arrests','Fontsize',16)

