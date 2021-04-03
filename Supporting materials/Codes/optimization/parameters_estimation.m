%script to optimize a functional of cost through CMA-ES

function [parameters_optim,cout] = parameters_estimation(varargin)

%INPUT:
%varargin: either empty or =
%support,repartition,data,repet,LoBounds,UpBounds

if nargin == 0 % called for a single estimation
    %% Preparation of data and initial distribution
    
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
    
    data(data == 0) = []; % in case the first long cycle recorded is the first cycle of the lineage
    
    data = flip(sort(data))';
    
    
    
    %% Monte-Carlo parameter
    
    repet = 250;
    
    
    %% Bounds, scaling and initial parameters
    
    LoBounds = zeros(2,1);
    UpBounds = zeros(2,1);
    
    i = 1;
    LoBounds(i) = 0.05;
    UpBounds(i) = 5;
    i = i+1;
    LoBounds(i) = 0.05;
    UpBounds(i) = 5;
    
    
else
    
    support = [varargin{1}];
    repartition = [varargin{2}];
    data = [varargin{3}];
    repet = [varargin{4}];
    LoBounds = [varargin{5}];
    UpBounds = [varargin{6}];
    
end


scaled_param0 = (log10(rand(size(LoBounds))) - log10(LoBounds))./(log10(UpBounds) - log10(LoBounds));



%% Preparation of CMAES
%%%%%%%%%%%%%minimisation%%%%%%%%%%%%%%%%%%%%
options.TolFun=1e-4;
options.TolX=1e-4;
%options.MaxFunEvals=20000;
options.MaxIter=20000;
options.LogPlot = 'off';

options.LBounds = zeros(size(LoBounds));
options.UBounds = ones(size(UpBounds));

sigma=[];

tic
[scaled_parameters_optim,cout,~,~,~,~]= cmaes('fitness_mean_error_exp',scaled_param0,sigma,options,repet,data,support,repartition,LoBounds,UpBounds);
toc

%% Plots

% figure;
plot(options.LBounds,1:numel(scaled_parameters_optim),'o',scaled_parameters_optim,1:numel(scaled_parameters_optim),'+',options.UBounds,1:numel(scaled_parameters_optim),'o','MarkerSize',12,'LineWidth',1.5)
hold on

%% Rescaling

parameters_optim = LoBounds.*((UpBounds./LoBounds).^scaled_parameters_optim);

end
