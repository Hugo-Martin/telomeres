function [fitted_curve,gof,robust,algo] = custom_best_fit(support,to_be_fitted,fitfun,x0,tol,maxiter,maxfunevals)

% INPUT: 
% support: x of data to fit, column vector
% to_be_fitted: y of data to fit, column vector
% fitfun: the custom equation, specified as a function handle
% x0: strating point
% tol: tolerance
% maxiter: maximum iterations
% maxfunevals: maximum function evaluations

% OUTPUT:
% fitted_curve: cfit objet that contains the information needed about the
% fitted curve
% gof: goodness of fit

type_of_fit = {'Off','Trust-Region';'Off','Levenberg-Marquardt';'LAR','Trust-Region';'LAR','Levenberg-Marquardt';'Bisquare','Trust-Region';'Bisquare','Levenberg-Marquardt'};
    
all_fits = cell(size(type_of_fit,1),3);
% 'Method','NonlinearLeastSquares',
for j = 1:size(type_of_fit,1)
    [f_c,gof_loc] = fit(support,to_be_fitted,fitfun,'StartPoint',x0,'Robust',type_of_fit{j,1},'Algorithm',type_of_fit{j,2},'TolFun',tol,'MaxIter',maxiter,'MaxFunEvals',maxfunevals);
    all_fits{j,2} = f_c;
    all_fits{j,3} = gof_loc;
    all_fits{j,1} = gof_loc.adjrsquare;
end

[~,i] = max([all_fits{:,1}]);

fitted_curve = all_fits{i,2};
gof = all_fits{i,3};
robust = type_of_fit{i,1};
algo = type_of_fit{i,2};




