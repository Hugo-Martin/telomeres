% Script to model the osnet of the first long cycle as a sequence of Bernouilli trials with variable probabilities p(j) at generation j

close all

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

fig_properties

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;

min_length_traj=3;

DOX=1; %1 to take only cycles after DOX addition, else all

[cell_cycles_lengths, indeces, lineages_length, ~]=extraction(data_exp,DOX,min_length_traj); % extraction of data


display_plots = 0; % to display plots for each value of the threshold

tried_threshhold = 13:30;

vect_p = zeros(size(tried_threshhold));


%% Loop on the value of the threshold

for i = 1:numel(tried_threshhold)
    
    threshold = tried_threshhold(i);
    
%% Preparation of data

data = [];
no_LC = [];

for j = 1:numel(data_exp)
    v = cell_cycles_lengths(indeces == j);
    if any(v > threshold)
        if any(v(find(v > threshold,1):end) <= threshold)
            data = [data, find(v > threshold,1)];
        end
    else
        no_LC = [no_LC, numel(v)]; % nb of generations in a lineage without long cycles
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

x0 = [50 10];

fitfun = fittype( @(a,b,x) 1./(1+exp(-(x-a)./b)));

tol = 1e-16;
maxiter = 5000;
maxfunevals = 8000;


[fitted_curve,gof,~,~] = custom_best_fit((1:numel(P_realistic))',P_realistic',fitfun,x0,tol,maxiter,maxfunevals);

c = coeffvalues(fitted_curve);
a = c(1); b = c(2);

fit_sigm = 1./(1 + exp(-((1:max(data)) - a)/b));


%% Creation of the probability mass function


pmf_sigm = fit_sigm(1);


for j = 2:max(data)
    pmf_sigm = [pmf_sigm, fit_sigm(j)*prod(1-fit_sigm(1:(j-1)))];
end



%% Statistical tests

[~, vect_p(i), ~,~,~] = sequence_bernouilli(data,pmf_sigm);


end

%% Plot

figure;
plot(tried_threshhold,vect_p,'--*',tried_threshhold,0.05*ones(1,numel(tried_threshhold)),'LineWidth',1.5,'MarkerSize',15)
xlabel('Value of the threshold')
ylabel('p-value')
legend({'Logistic fit'})