function [seq,type] = structure(div, threshold)

% function that identify the 'structure' of a lineage: number of cycles in
% sequences in each state (short or long)

%OUTPUT
%seq: vector made of the lengths of the different sequences (of short or long cycles), in
%generations, in each lineages

%type : vector of the same size than seq made of 0 (for sequences of short
%cycles) or 1 (long cycles)

%INPUT
%div : vector representing a lineage, each entry is the duration of a cycle
%threshold : difference short/long

bool = div>threshold;
state = bool(1); %current state

seq = [];
type = state;

old = 1; %counter for the sequences
new = 1;



while isempty(new) == 0 %while the end of the cycle is reached
    new = find(bool(old:end) ~= state,1); %the first index of a cycle with a different state?
    if isempty(new) == 0
        seq = [seq; new-1];
        old = old + new -1;
        state = bool(old);
        type = [type; state];
    end
end

seq = [seq; numel(div)-old+1]; %last update of seq


