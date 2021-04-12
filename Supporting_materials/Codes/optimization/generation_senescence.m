% function to generate a generation of senescence for a lineage

function [nn,len] = generation_senescence(a,b,support,repartition)
% [nn,len]
% INPUT:
% a,b: parameters
% support: the support of the initial distribution
% repartition: its cumulative distribution function

% OUTPUT:
% nn: generation of senescence for one lineage
% len: length of the shortest telomere at signaling


n = 32; % nb of telomeres

L = random_generator(n,support,repartition);


L = round(L*2)/2; % rounded at the closest 0.5, as in support

len0=min(L,[],'all');
k=1;
            
            L = reshape(L(randperm(numel(L))),2,n/2); % random_generation_initial_length provides a sorted list of values
            test = rand;

           while and(((b*exp(-a*min(min(L)))) <= test), all(L > 0)) % while none of the telomeres are short enough to signal senescence
               len0 = min(min((L)));
               B = unidrnd(2,1,n/2)-1;% selection of which telomere in each pair will be shortened
                R = 4*ones(1,n/2) + random('Discrete Uniform',6,1,n/2);% the length which is cut from telomeres
                L = L - [R.*B; R.*(~B)];
                k = k+1;    
                test = rand;
            end
            

nn=k;
len=len0;
end
