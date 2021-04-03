function [h, p, st,edges,binned_expected] = sequence_bernouilli(data,pmf)


% Function that first bins the data and expected counts so that at most one class
% of expected counts has fewer than 5 counts and second performs a chi
% square goodness of fit test

% OUTPUT:
%
% h: 1 rejection, 0 no rejection (standard 0.05 threshold)
% p: p value
% st: statistics of the test

% INPUT:
%
% data: data processed as in FLC_Bernouilli
% pmf: associated probability mass function, once the 'instantaneous'
% probability of forst long cycle is fitted on a S curve

edges = 0.5;
expected = numel(data)*pmf;
local_bin_value = 0;
local_counter = 0;
binned_expected = [];

while ~isempty(expected) % while expected is not empty
    local_bin_value = local_bin_value + expected(1);
    local_counter = local_counter+1;
    if local_bin_value >= 5 % usual minimum value for (most of) bins
        edges = [edges, edges(end) + local_counter];
        binned_expected = [binned_expected, local_bin_value];
        local_bin_value = 0;
        local_counter = 0;
    end
    expected(1) = [];
end

if local_counter ~= 0 % if the last entry in the while loop did not "end" a bin
    edges = [edges, numel(pmf)+0.5];
    binned_expected = [binned_expected, local_bin_value];
end

[h,p,st] = chi2gof(data,'Edges',edges,'Expected',binned_expected,'NParams',2); % two parameters are determined from the data in the S curve
