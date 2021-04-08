%script to assess the variations in cell cycle duration in wild type
%lineages



close all

min_length_traj=3;

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

fig_properties

DOX=0; %1 to take only cycles after DOX addition, else all


trcritere31=load('TelomerasePositive1.mat'); %31 lineages
data_exp=trcritere31.OrdtryT528sansdox;

[cell_cycles_lengths, indeces, lineages_length, ~]=extraction(data_exp,DOX,min_length_traj); % extraction of data


data = zeros(numel(data_exp),max(lineages_length)-3); % table to store the cell cycle lengths

for j = 1:numel(data_exp)
    v = cell_cycles_lengths(indeces == j); % lineage j
    v = v(2:end-2); % censoring of the first and the two last cycles
    data(j,1:numel(v)) = v'; % line j is filed with the lengths
end


quant = [0.25 0.5 0.75]; % the quartiles

Q = zeros(numel(quant)+1, size(data,2));

for j = 1:size(data,2)
    col = data(:,j);
    col = col(col ~= 0); % exclusion of the zeros in each column
    Q(:,j) = [quantile(col,quant)'; numel(col)];
end



figure
subplot(2,1,1)
color(1,:)=[0.7 0.7 0.7];
color(3,:)=[0.7 0.7 0.7];
color(2,:)=[0 0 0];

for j = 1:numel(quant)
    plot(1:size(data,2),10*Q(j,:),'color',color(j,:),'LineWidth',3,'MarkerSize',30)
    hold on
end

ylim([0 400])

ax = gca;
ax.FontSize = 22;
ax.LineWidth=3;
ax.TickDir = 'out';

xlabel(label_gene)
ylabel(label_celldur)

legend({'1st quartile','Median','3rd quartile'},'FontSize',22,'location','northwest')
title('Evolution of cell cycle durations for the WT as generations pass by')

hold off


generations = [5 25 45];


limite = 30;
subplot(2,1,2)


for j = 1:numel(generations)
    v = data(data(:,generations(j))~=0,generations(j)); % columns = generations 5, 25 and 45, excluding the values 0
    v = [zeros(1,min(v)-1), histcounts(v,'BinMethod','integers')]; % adding 0 at the begining of each rows, corresponding to non existant durations for the considered generation
    v = [v, zeros(1,limite-numel(v))];
    v = v/sum(v);
    v = cumsum(v);
    plot(10*(1:limite),v(1:limite),'LineWidth',3)
    hold on
end


ax = gca;
ax.FontSize = 22;
ax.LineWidth=3;
ax.TickDir = 'out';

xlabel(label_celldur)
ylabel('Cumulative frequency')
legend({'Generation 5','Generation 25','Generation 45'},'FontSize',22,'location','southeast')
%title({'Comparison of the cumulative frenquency of the durations of the cell cycles', 'for three generations of wild type yeasts'},'FontSize',26)

hold off
set(gcf, 'Color', [1,1,1]);


%savefig('Fig1.fig') % uncomment if you want to save the figure