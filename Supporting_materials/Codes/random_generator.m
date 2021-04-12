% Function that generate a sample of realizations with a given cumulative
% distribution function.

function mysample = random_generator(n,support,repartition)

% INPUT:
% n: number of initial telomere length to generate
% support: support of the distribution of initial lengths, as a column
% vector
% repartition: cumulative mass function of initial lengths, as a
% column vector /!\ need to be sorted /!\


% OUPUT:
% mysample: sorted n-sample of realizations of the provided law,
% as a column vector



random_sample = sort(unifrnd(0,1,n,1));

Large = random_sample(random_sample>max(repartition));

if ~isempty(Large) %If there are values larger than the maximum of the CDF
    random_sample = setdiff(random_sample,Large); %the loop is runned on the number of small enough values
end

Small = random_sample(random_sample<min(repartition));

if ~isempty(Small) %Similar than before, but with smaller values
    random_sample = setdiff(random_sample,Small);
end


if isempty(random_sample) % if there are only 'small' and 'large' values
    mysample = [support(1)*ones(numel(Small),1);support(end)*ones(numel(Large),1)];
else
    mysample = zeros(size(random_sample));
    mysample(1) = find(repartition>=random_sample(1),1);
    
    for c = 2:numel(random_sample)
        mysample(c) = find(repartition >= random_sample(c),1); 
    end
    mysample = [support(1)*ones(numel(Small),1);support(mysample);support(end)*ones(numel(Large),1)];
end


