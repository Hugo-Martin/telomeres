% Script to model the osnet of the first long cycle as a sequence of Bernouilli trials with variable probabilities p(j) at generation j

close all
fig_properties

addpath('./Supporting material/Codes')
addpath('./Supporting material/Data')

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;

min_length_traj=3;

DOX=1; %1 to take only cycles after DOX addition, else all

[cell_cycles_lengths, indeces, lineages_length, ~]=extraction(data_exp,DOX,min_length_traj); % extraction of data


threshold = 18;


%% Preparation of data

data = [];
no_LC = [];

for j = 1:numel(data_exp)
    v = cell_cycles_lengths(indeces == j);
    if any(v > threshold) % if there is a long cycle
        if any(v(find(v > threshold,1):end) <= threshold) % if any short cycle after the first long
            data = [data, find(v > threshold,1)];
        end
    else
        no_LC = [no_LC, numel(v)]; % durée en générations des lignées sans cycle long
    end
end

%for each subset, computation of the vector that contains the total number
%of lineages that do not have had its first long cycle at each generation


concat = [data, no_LC];

remaining_pool = zeros(1,max(data));

for k = 1:numel(data)
    remaining_pool = remaining_pool + [ones(1,concat(k)),zeros(1,numel(remaining_pool) - concat(k))];
end

realisations = zeros(1,max(concat));

for k = 1:numel(data)
    realisations(data(k)) = realisations(data(k)) + 1;
end


P = realisations./remaining_pool;



min_lignees = 10; % minimum number of data to compute the 'instantaneous' probability

P_realistic = P(remaining_pool >= min_lignees);


%% Fitting on a S curve

%
x0 = [50 10];

fitfun = fittype( @(a,b,x) 1./(1+exp(-(x-a)./b)));

tol = 1e-16;
maxiter = 5000;
maxfunevals = 8000;


[fitted_curve,gof,~,~] = custom_best_fit((1:numel(P_realiste_internal))',P_realiste_internal',fitfun,x0,tol,maxiter,maxfunevals);

c = coeffvalues(fitted_curve);
a_internal = c(1); b_internal = c(2);
%a = 59.07; b = 14.3;
fit_sigm = 1./(1 + exp(-((1:max(data)) - a)/b));


%% Creation of the probability mass function


pmf_sigm = fit_sigm(1);

for j = 2:max(data)
    pmf_sigm = [pmf_sigm, fit_sigm(j)*prod(1-fit_sigm(1:(j-1)))];
end



%% Plots

repet = 1000;
simulated_data = zeros(numel(data),repet);
for j = 1:repet
    mysample = random_generator(numel(data),(1:numel(pmf_sigm))',cumsum(pmf_sigm'));
    simulated_data(:,j) = mysample;
end

quant = quantile(simulated_data,[0.025, 0.975],2);

figure(1);

hold on

quant_2=zeros(numel(data)+1,2);
quant_2(1:numel(data),:)=quant;
quant_2(numel(data)+1,:)=max(data);
ar1=area(quant_2(:,1),1:numel(data)+1,'FaceColor',[191/255 220/255 234/255],'HandleVisibility','off');
ar2=area(quant_2(:,2),1:numel(data)+1,'FaceColor','w','HandleVisibility','off');
ar1.EdgeColor='white';
ar2.EdgeColor='white';

h2=plot(quant(:,1),1:numel(data),'Color',[100/255 100/255 200/255],'LineWidth',3);
h3=plot(quant(:,2),1:numel(data),'Color',[100/255 100/255 200/255],'LineWidth',3);

%116/255 208/255 241/255
%h3=plot(quant(:,2),1:numel(data),'Color',[82/255 145/255 225/255],'LineWidth',3);
%title('Data VS the series of Bernoulli model')
h1=plot(sort(data),1:numel(data),'k','LineWidth',3);
ax_properties
xlim([0 max(data)]);

xlabel('Generation of the first arrest')
ylabel('Index of lineages')

legend({'2.5% quantile','97.5% quantile','Data'},'Location','Northwest')
%savefig('Fig2A.fig')


figure(2)
plot(1:numel(P_realistic),P_realistic,'o','MarkerSize',15,'MarkerFaceColor',ColorData2,'MarkerEdgeColor',ColorData2)
hold on
plot(1:numel(fit_sigm),fit_sigm,'LineWidth',3,'Color',ColorFit1)
xlabel('Generation')
ylabel('Probability')
legend('Data','Logistic fit','location','northwest')

ax_properties
%savefig('Fig2B.fig')
%title('Probability of occurence of the first long cycle at each generation')


phat = mle(data,'distribution','geo');

figure(3);

histogram(data,'BinMethod','integers')%,'Normalization','probability');
hold on
plot(1:numel(pmf_sigm),numel(data)*pmf_sigm,'o',1:numel(pmf_sigm),numel(data)*pdf('Geometric',1:numel(pmf_sigm),phat),'o')
hold off
legend('Data','Logistic parameter','Constant parameter')
xlabel('Generation')
title('Generation of occurence of the first internal long cycle, internal')
%savefig('Fig2C.fig')


