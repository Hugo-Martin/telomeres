% Function qui tire au hasard numel(x) entiers
% selon la fonction de répartition donnée. Utilise crucialement que le vecteur
% repartition est trié.


function mysample = random_generator(n,support,repartition)

% INPUT:
% n: number of initial telomere length to generate
% support: support of the distribution of initial lengths, as a column
% vector
% repartition: cumulative mass function of initial lengths, as a
% column vector


% OUPUT:
% mysample: sorted n-sample of realizations of the provided law,
% as a column vector



random_sample = sort(unifrnd(0,1,n,1));

grandes = random_sample(random_sample>max(repartition));

if ~isempty(grandes) %S'il y a des valeurs plus grandes que le maximum de repartition
   % disp([num2str(numel(grandes)), ' valeurs plus grandes que le maximum de repartition'])
    random_sample = setdiff(random_sample,grandes); %on va lancer la boucle sur les valeurs suffisamment petites
end

petites = random_sample(random_sample<min(repartition));

if ~isempty(petites) %S'il y a des valeurs plus grandes que le maximum de repartition
  %  disp([num2str(numel(petites)), ' valeurs plus petites que le minimum de repartition'])
    random_sample = setdiff(random_sample,petites); %on va lancer la boucle sur les valeurs suffisamment grandes
end


if isempty(random_sample) % if there are only 'small' and 'large' values
    mysample = [support(1)*ones(numel(petites),1);support(end)*ones(numel(grandes),1)]; % On ajoute autant de valeurs extrêmes que de valeurs trop petites et trop grandes
else
    mysample = zeros(size(random_sample));
    mysample(1) = find(repartition>=random_sample(1),1);
    
    for c = 2:numel(random_sample)
        mysample(c) = find(repartition >= random_sample(c),1); 
    end
    mysample = [support(1)*ones(numel(petites),1);support(mysample);support(end)*ones(numel(grandes),1)]; % On ajoute autant de valeurs extrêmes que de valeurs trop petites et trop grandes
end


