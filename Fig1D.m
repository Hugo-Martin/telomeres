close all

addpath('./Supporting_materials/Codes')
addpath('./Supporting_materials/Data')

min_length_traj=3;

fig_properties

tr=load('TelomeraseNegative.mat');
data_exp=tr.OrdtryT528total160831;


DOX_exp=1; %1 si on veut après l'addition de DOX, 0 si on prend tout.

[cycles_duration, indeces, ~, ended]=extraction(data_exp,DOX_exp,min_length_traj);

internal = [];
senescence = [];

threshold = 18;

for j = 1:numel(data_exp)
    v = cycles_duration(indeces == j);
    if (ended(j) == 1) && (v(end) > threshold) && (any(v <= threshold))%= senescent lineage that contains at least one short cycle
        senescence = [senescence, numel(v) - find(flip(v) <= threshold,1) + 2];
    end
    f = find(v > threshold,1);
    if ~isempty(f)
        if any(v(f:end) >= threshold)
            internal = [internal, f];
        end
    end
end


figure(1)
[F_sene,X_sene]=ecdf(senescence);
[F_intern,X_intern]=ecdf(internal);
plot(X_intern,F_intern,X_sene,F_sene,'LineWidth',linewidth)
ax_properties
xlabel('Generation of the 1st Long Cycle')
ylabel('Empirical cumulated distribution function')
%savefig('Fig1C_TTMarch21.fig')