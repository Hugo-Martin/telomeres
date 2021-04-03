% Script to visualize the simulations and data

%close all

format long

addpath('./Supporting material/Codes')
addpath('./Supporting material/Data')
addpath('./Supporting material/Codes/optimization')

fig_properties
%% Preparation of data and initial distribution


load('etat_asymp_val_juillet');
support = etat_asymp_val_juillet(1:numel(etat_asymp_val_juillet)/2); % support
repartition = cumsum(etat_asymp_val_juillet(numel(etat_asymp_val_juillet)/2+1:end)); % cumulative mass function

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;

[dureesdiv,indice,longueur_lignees, lignees_terminees]=extraction(data_exp,1,3);

% First arrest, 187 lineages
threshold = 18; % seuil long/court

data_first_arrest = [];


%Arrests that are not the begining of senescence only

data_first_arrest_non_senescent = [];


for j = 1:numel(data_exp)
    v = dureesdiv(indice == j);
    first_long = find(v>threshold,1);
    if (~isempty(first_long)) && (any(v(first_long:end) <= threshold))
        data_first_arrest_non_senescent=[data_first_arrest_non_senescent, first_long];
    end
end

%data_first_arrest_non_senescent(data_first_arrest_non_senescent == 0) = []; % in case the first long is the first cycle recorded

data_first_arrest = flip(sort(data_first_arrest))';
data_first_arrest_non_senescent = flip(sort(data_first_arrest_non_senescent))';

%% Parameters for p(L)=b*exp(-a*L)

%b=0.17;
%a=0.0173;

a = 0.023;
b = 0.276;

%% Simulations

Nb_simu = 1000;



simu_first_arrest = zeros(1,Nb_simu*numel(data_first_arrest));
simu_first_arrest_non_senescent = zeros(1,Nb_simu*numel(data_first_arrest_non_senescent));

len=zeros(1,Nb_simu*numel(data_first_arrest));
len_ns=zeros(1,Nb_simu*numel(data_first_arrest_non_senescent));
vect_len_first_arrest = simu_first_arrest;
vect_pos_first_arrest = simu_first_arrest;
vect_len_first_arrest_ns = simu_first_arrest_non_senescent;
vect_pos_first_arrest_ns = simu_first_arrest_non_senescent;

% parfor k = 1:numel(simu_first_arrest)
%     [simu_first_arrest(k) len(k)]= generation_senescence_expo(a,b,support,repartition);
% end

parfor k = 1:numel(simu_first_arrest_non_senescent)
    [simu_first_arrest_non_senescent(k), len_ns(k)]= generation_senescence_expo(a,b,support,repartition);
end

%% Plot


%simulated_data_first_arrest = sort(reshape(simu_first_arrest,numel(data_first_arrest),Nb_simu));
%quant_first_arrest = quantile(simulated_data_first_arrest,[0.025, 0.975],2);

simulated_data_first_arrest_ns = sort(reshape(simu_first_arrest_non_senescent,numel(data_first_arrest_non_senescent),Nb_simu));
quant_first_arrest_ns = quantile(simulated_data_first_arrest_ns,[0.025, 0.975],2);

% figure;
% %subplot(1,3,1)
% hold on
% %plot(flip(data_first_arrest),1:numel(data_first_arrest),quant_first_arrest(:,1),1:numel(data_first_arrest),'-*',quant_first_arrest(:,2),1:numel(data_first_arrest),'-*',min(simulated_data_first_arrest,[],2),1:numel(data_first_arrest),'-o',max(simulated_data_first_arrest,[],2),1:numel(data_first_arrest),'-o')
% %title('First arrest')
% hold on
% quant_2=zeros(numel(data_first_arrest)+1,2);
% quant_2(1:numel(data_first_arrest),:)=quantile(simulated_data_first_arrest,[0.05,0.95],2);
% quant_2(numel(data_first_arrest)+1,:)=max(data_first_arrest);
% ar1=area(quant_2(:,1),1:numel(data_first_arrest)+1,'FaceColor',[191/255 220/255 234/255],'HandleVisibility','off')
% ar2=area(quant_2(:,2),1:numel(data_first_arrest)+1,'FaceColor','w','HandleVisibility','off');
% ar1.EdgeColor='white';
% ar2.EdgeColor='white';
% 
% h2=plot(quant_2(1:numel(data_first_arrest),1),1:numel(data_first_arrest),'Color',[100/255 100/255 200/255],'LineWidth',3);
% h3=plot(quant_2(1:numel(data_first_arrest),2),1:numel(data_first_arrest),'Color',[100/255 100/255 200/255],'LineWidth',3);
% 
% %116/255 208/255 241/255
% %h3=plot(quant(:,2),1:numel(data),'Color',[82/255 145/255 225/255],'LineWidth',3);
% %title('Data VS the series of Bernoulli model')
% h1=plot(sort(data_first_arrest),1:numel(data_first_arrest),'k','LineWidth',3);
% ax_properties
% xlim([0 max(data_first_arrest)]);
% 
% xlabel('Generation of the first arrest')
% ylabel('Index of lineages')
% 
%  legend({'5% quantile','95% quantile','Data'},'Location','Northwest')
% savefig('Fig17_exp.fig')
% 

x = 0:0.1:250;
y_first_arrest = b.*exp(-a*x);



figure;
plot(x,y_first_arrest,'LineWidth',linewidth)
%title('Instantaneous probability of signalling')
ax_properties
xlabel('Telomere length')
ylabel('Probability of signalling senescence')
savefig('Fig18_exp.fig')

% figure;
% histogram(len,'BinMethod','integers','normalization','probability')
% xlabel('Length of the signalling telomere')
% ylabel('Frequency')
% ax_properties
% %title('Probability distribution of the telomere length at signaling')
% savefig('Fig16_exp.fig')
% 

%%%%FIGURES BIS: SEULEMENT POUR LES 1ers ARRETS NON SENESCENTS

figure;
%subplot(1,3,1)
hold on
%plot(flip(data_first_arrest),1:numel(data_first_arrest),quant_first_arrest(:,1),1:numel(data_first_arrest),'-*',quant_first_arrest(:,2),1:numel(data_first_arrest),'-*',min(simulated_data_first_arrest,[],2),1:numel(data_first_arrest),'-o',max(simulated_data_first_arrest,[],2),1:numel(data_first_arrest),'-o')
%title('First arrest')
hold on
quant_2_ns=zeros(numel(data_first_arrest_non_senescent)+1,2);
quant_2_ns(1:numel(data_first_arrest_non_senescent),:)=quantile(simulated_data_first_arrest_ns,[0.05,0.95],2);
quant_2_ns(numel(data_first_arrest_non_senescent)+1,:)=max(data_first_arrest_non_senescent);
ar1=area(quant_2_ns(:,1),1:numel(data_first_arrest_non_senescent)+1,'FaceColor',[191/255 220/255 234/255],'HandleVisibility','off');
ar2=area(quant_2_ns(:,2),1:numel(data_first_arrest_non_senescent)+1,'FaceColor','w','HandleVisibility','off');
ar1.EdgeColor='white';
ar2.EdgeColor='white';

h2=plot(quant_2_ns(1:numel(data_first_arrest_non_senescent),1),1:numel(data_first_arrest_non_senescent),'Color',[100/255 100/255 200/255],'LineWidth',3);
h3=plot(quant_2_ns(1:numel(data_first_arrest_non_senescent),2),1:numel(data_first_arrest_non_senescent),'Color',[100/255 100/255 200/255],'LineWidth',3);

%116/255 208/255 241/255
%h3=plot(quant(:,2),1:numel(data),'Color',[82/255 145/255 225/255],'LineWidth',3);
%title('Data VS the series of Bernoulli model')
h1=plot(sort(data_first_arrest_non_senescent),1:numel(data_first_arrest_non_senescent),'k','LineWidth',3);
ax_properties
xlim([0 max(data_first_arrest_non_senescent)]);

xlabel('Generation of the first arrest')
ylabel('Index of lineages')

 legend({'2.5% quantile','97.5% quantile','Data'},'Location','Northwest')
savefig('Fig17_expBIS.fig')




figure;
histogram(len_ns,'BinMethod','integers','normalization','probability')
xlabel('Length of the signalling telomere')
ylabel('Frequency')
ax_properties
%title('Probability distribution of the telomere length at signaling')
savefig('Fig16_expBIS.fig')
