% function to generate a generation of senescence for a lineage

function [nn,len] = generation_senescence_expo(a,b,support,repartition)
% [nn,len]
% INPUT:
% a,b: parameters
% support: the support of the initial distribution
% repartition: its cumulative distribution function

% OUTPUT:
% nn: generation of senescence for one lineage
% len: length of the shortest telomere at signaling


n = 32; % nb of telomeres

L = random_generation_initial_length(n,support,repartition);


L = round(L*2)/2; % rounded at the closest 0.5, as in support

len0=min(L,[],'all');
k=1;
            
            L = reshape(L(randperm(numel(L))),2,n/2); % random_generation_initial_length provides a sorted list of values
            test = rand;

           while and(((b*exp(-a*min(min(L)))) <= test), all(L > 0)) % while none of the telomeres are short enough to signal senescence
               len0 = min(min((L)));
               B = unidrnd(2,1,n/2)-1;% le télomère de chaque pair qui sera réduit, vecteur de 0 et 1
                R = 4*ones(1,n/2) + random('Discrete Uniform',6,1,n/2);% la quantité dont seront réduits ces télomères
                L = L - [R.*B; R.*(~B)];
                k = k+1;    
                test = rand;
            end
            

nn=k;
len=len0;
end
