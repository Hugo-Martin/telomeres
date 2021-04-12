%Study of the onset of the second long cycle

close all

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

min_length_traj=3;

fig_properties

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;


tested_threshold = 9:40;

DOX_exp=1;

[cycles_duration, indeces, ~, ~]=extraction(data_exp,DOX_exp,min_length_traj);


pval = zeros(size(tested_threshold));

proportion_lineages_with_two_consecutive_long_cycles = zeros(size(tested_threshold));
prop_one_div_intern=zeros(size(tested_threshold));

for i = 1:numel(tested_threshold)
    threshold = tested_threshold(i);

data = [];

for j = 1:numel(data_exp)
    v = cycles_duration(indeces == j);
    f = find(v > threshold, 2);
    if (numel(f) == 2) && (any(v(f(2):end) <= threshold))
        data = [data,f];
    end
end

proportion_lineages_with_two_consecutive_long_cycles(i) = sum(data(2,:) - data(1,:) == 1)./size(data,2);



[~, ~, ~, pval(i)] = test_geom(data(2,:)-data(1,:), 1, 0, size(data,1),0);


if threshold == 18 % to illustrate
    
ax_properties
ylim([0 max(data(1,:))])
xlim([0 max(max(data(1,:)))])

savefig('Fig11.fig')

figure(1);
plot(data(1,:),data(2,:)-data(1,:),'o','MarkerSize',markersize,'MarkerFaceColor',ColorData2,'MarkerEdgeColor',ColorData2)

xlabel('Generation of the 1st arrest')
ylabel('Distance between the 1st and 2nd arrests')


ax_properties
ylim([0 10*ceil(max(data(2,:) - data(1,:))/10)])
xlim([0 10*ceil(max(data(1,:))/10)])

savefig('Fig4B_TTMarch21.fig')


figure(2)
histogram(data(2,:)-data(1,:),'BinMethod','integers','EdgeColor',ColorData1,'FaceColor',ColorData2)
xlabel('Number of generations between the first two long cycles')
ylabel('Count')
ax_properties
savefig('Fig4A_TTMarch21.fig')
end

end
% 
figure(3);
plot(tested_threshold*10,pval,'o','MarkerSize',15,'MarkerFaceColor',ColorFit1,'MarkerEdgeColor',ColorFit1)
hold on
plot(tested_threshold*10,0.05*ones(size(pval)),'LineWidth',3,'Color',ColorFit1)
xlabel('Threshold (min)')
ylabel('p-value')
ax_properties
ylim([0 1])
xlim([min(tested_threshold)*10 max(tested_threshold)*10])
savefig('Fig10.fig')


figure(4)
plot(tested_threshold*10,proportion_lineages_with_two_consecutive_long_cycles,'o','MarkerSize',markersize,'MarkerFaceColor',ColorData2,'MarkerEdgeColor',ColorData2)
xlabel('Threshold (min)')
ylim([0 1])
xlim([min(tested_threshold*10) max(tested_threshold*10)])
ylabel('Percentage of two consecutive first long cycles')
ax_properties
savefig('Fig12.fig')

