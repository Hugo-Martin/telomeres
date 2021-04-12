function [phat, pci, h, pval, st] = test_geom(values, expected_min, plotting, nbreal,Emin)

% Function that tests if a given distribution is consistent with a
% geometric distribution

% If there are too few categories (e.g. chi-squared test, less than three categories), the displayed value is 0. 

% OUTPUT:
%
% phat : estimated parameter with maximum likelyhood estimate
% pci : 95% confidence interval for the estimated parameter
% pval: p-value given by the test

% INPUT:
%
% values : list of values
% expected_min : 0 or 1, depending on the minimum expected value for the
% studied dataset
% plotting : 1 to display a graphic display, 0 otherwise



if expected_min == 1
    shift = 1;
elseif expected_min == 0
    shift = 0;
else
    disp('Wrong minimal expected value')
end

% values <- values - shift : in MATLAB, a géometric law can have the value 0

[phat, pci] = mle(values - shift,'distribution','geo'); % maximum likelihood estimate

histcnts = histcounts(values - shift,'BinMethod','integers');


pd = makedist('NegativeBinomial','R',1,'p',phat); % making of a geometric distribution with the estimated parameter

bins = 0:(numel(histcnts)-1);

expected = sum(histcnts) * pdf(pd,bins);

[h,pval,st] = chi2gof(bins,'Ctrs',bins,'Frequency',histcnts,'Expected',expected,'NParams',1,'EMin',Emin);


    
%% Graphic representation

    if plotting == 1 
        
        M = max(values+shift);
        
        table = zeros(nbreal,M+1-shift);
        
        for k = 1:nbreal
            rand = geornd(phat,1,numel(values))+1; % +1 because in matlab a geometric law starts at 0
            r = histcounts(rand);
            if numel(r)<M+1-shift
                r = padarray(r',M+1-shift-numel(r),'post')'; % adaptation of the size of h
                table(k,:) = r;
            elseif numel(r)>M+1-shift
            else
                table(k,:) = r;
            end
        end
        
        Q = quantile(table,[0.025 0.975],1);
        
        
        figure;
        histogram(values,'BinMethod','integers');
        hold on
        plotting(shift:size(Q,2)-1+shift,Q(2,:),'*',shift:size(Q,2)-1+shift,Q(1,:),'*')
        xlabel('Realizations','FontSize',18)
        ylabel('Counts')
        legend('Experimental data','Lower bound of the 95% quantile','Upper bound of the 95% quantile')
        hold off
        
    end