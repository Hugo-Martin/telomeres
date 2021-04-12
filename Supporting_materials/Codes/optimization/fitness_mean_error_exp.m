%Compute the error in norm 2 for a given couple of parameters (mu,sigma).

function e = fitness_mean_error_exp(param,repet,sorted_data,support,repartition,varargin)
 
% INPUT:
%param: vector of parameters [a,b]
% repet: number of repetitions of the experiment for the Monte-Carlo mean
% sorted_data: sorted list of experimental data as a sorted column vector,
% with smaller values at then end
% support: support of the initial distribution
% repartition: cumulative mass function of the initial distribution
% varargin: either empty or = LoBounds,UpBounds


% OUTPUT:
% e: error associated to a couple of parameters (a,b)

%disp(nargin)

if nargin == 7 % called with scaled data
    LoBounds = [varargin{1}];
    UpBounds = [varargin{2}];
    a = LoBounds(1)*((UpBounds(1)/LoBounds(1))^param(1));
    b = LoBounds(2)*((UpBounds(2)/LoBounds(2))^param(2));
elseif nargin == 5 % called with unscaled data
    a = param(1);
    b = param(2);
else
    disp('Wrong number of input arguments, must be either 4 or 6')
end


nn=zeros(numel(sorted_data)*repet,1);


    parfor j=1:numel(sorted_data)*repet
        nn(j) = generation_senescence(a,b,support,repartition);
    end

nn=flip(sort(reshape(nn,numel(sorted_data),repet)));


repetde = repmat(sorted_data,1,repet);

e = sum(sum((nn-repetde).^2));


end