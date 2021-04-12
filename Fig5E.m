% Boxplots to study how long the arrests are


close all

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

fig_properties

min_length_traj=3;

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;



DOX_exp=1;



threshold = 18;

[cycle_lengths, indeces, ~, ended]=extraction(data_exp,DOX_exp,min_length_traj);

last_of_sene = [];
last_of_others = [];

others_of_sene = [];
others_of_others = [];


for j = 1:numel(data_exp)
    v = cycle_lengths(indeces == j);
    %% If the last cycle is long, the lineage is treated differently, depending on if the lineage is dead at the end of the experiment
    if v(end) > threshold % if the last cycle is long
        if ended(j) == 1 % if additionaly, the lineage is ended at the end of the experiment, thus senescent
            last_of_sene = [last_of_sene, v(end)];
            v(end) = []; % removing of the last cycle
            if v(end) > threshold % if the 'new last cycle' is long
                vect = flip(v);
                vect = vect(1:find(vect <= threshold, 1)-1);
                others_of_sene = [others_of_sene, vect'];
                v = v(1:end-numel(vect)); % removing of the remaining of the senescent phase
            end
        else % if the lineage is not ended at the end of the experiment but the last cycle is long, we cannot know for sure if it would have been the senescent phase if the experiment would have lasted longer
            v(numel(v) - find(flip(v) <= threshold, 1) + 2:end) = [];
        end
    end
    %% If the last cycle is short
    while isempty(v) == 0 % while all the lineage is not studied
        i = find(v > threshold,1);
        k = find(v(i:end) <= threshold,1);
        if  isempty(k) % if there is not short cycle after the i-th cycle (which is long)
            if ~isempty(i) % if there remain long cycles
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


boxplot([last_of_others, last_of_sene],[repmat({'Internal sequences'},1,numel(last_of_others)),repmat({'Senescence'},1,numel(last_of_sene))],'LabelVerbosity','minor')
title('Last cycle in the sequences of arrests','Fontsize',16)

figure(1920)
h1=boxplot([others_of_sene, last_of_sene, others_of_others,last_of_others],[repmat({'Previous cycles'},1,numel(others_of_sene)),repmat({'Last cycles'},1,numel(last_of_sene)),repmat({'Previous cycles'},1,numel(others_of_others)),repmat({'Last cycles'},1,numel(last_of_others))],'LabelVerbosity','minor');
set(h1(:,:),'LineWidth',linewidth)
set(h1(:,:),'Color',ColorData2)
ax_properties
ylim([0 300])
%savefig('Fig1920.fig')

figure(19)
h=boxplot([others_of_others, last_of_others],[repmat({'Previous cycles'},1,numel(others_of_others)),repmat({'Last cycles'},1,numel(last_of_others))],'LabelVerbosity','minor');
set(h(:,:),'LineWidth',linewidth)
set(h(:,:),'Color',ColorData2)
ax_properties
ylim([0 300])
title('Internal sequences of arrests')
%savefig('Fig19.fig')

figure
boxplot([others_of_others, others_of_sene],[repmat({'Other sequences'},1,numel(others_of_others)),repmat({'Senescence'},1,numel(others_of_sene))],'LabelVerbosity','minor')
title('Cycles before the last of sequences of arrests','Fontsize',16)

