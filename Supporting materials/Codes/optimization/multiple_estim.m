% script to estimate the parameters multiple times and store the results

%% Preparation of data and parameters

addpath('../')
addpath('../../Data')


load('etat_asymp_val_juillet');
support = etat_asymp_val_juillet(1:numel(etat_asymp_val_juillet)/2); % support
repartition = cumsum(etat_asymp_val_juillet(numel(etat_asymp_val_juillet)/2+1:end)); % cumulative mass function


tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;

[duration,indeces,~, ~]=extraction(data_exp,1,3);

% First arrest, non terminal arrests

threshold = 18;

data = [];

for j = 1:numel(data_exp)
    v = duration(indeces == j);
    if any(v > threshold)
        l = find(v >threshold,1);
        if any(v(l:end) < threshold)
            data = [data, l];
        end
    end
end


data = flip(sort(data))';



% Monte-Carlo parameter

repet = 250;


% Bounds, scaling and initial parameters

LoBounds = zeros(2,1);
UpBounds = zeros(2,1);

i = 1;
LoBounds(i) = 0.001;
UpBounds(i) = 1;
i = i+1;
LoBounds(i) = 0.01;
UpBounds(i) = 2;

%% Estmations

nb_estim = 10;

estimated = zeros(2,nb_estim);
cost = zeros(1,nb_estim);

fid = fopen('optim_exp.txt','w');
fprintf(fid,'%6s %6s %12s\r\n','a','b','fitness');
fclose(fid);

for i = 1:nb_estim
    [param,fitness] = parameters_estimation(support,repartition,data,repet,LoBounds,UpBounds);
    estimated(:,i) = param;
    cost(i) = fitness;
    fid = fopen('optim_exp.txt','a');
    fprintf(fid,'%2.6f %2.6f %2.6f\n',[param',fitness]);
end




%% Data analysis (identifiability)



fileID = fopen('optim_exp.txt','r');
M = fscanf(fileID,'%f %f %f',Inf);
M = reshape(M,[3,numel(M)/3]);
M = M';

min_fitness = min(M(:,3));

small_fitness = M(:,3) < 1.05*min_fitness; % the 5% smallest fitnesses

names = {'a';'b';'fitness'};

figure;
for j = 1:3
    arguments = M(small_fitness,j);
    subplot(1,3,j)
    if j < 3
        [f,xi] = ksdensity(arguments);
        plot(arguments,zeros(size(arguments)),'+',xi,f)
        xlim([LoBounds(j) UpBounds(j)]);
    else
        [f,xi] = ksdensity(arguments);
        plot(arguments,zeros(size(arguments)),'+',xi,f)
        xlim([min(arguments) max(arguments)]);
    end
    title(names{j})
end
    


