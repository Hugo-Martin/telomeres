function [cycles_length,indeces,ligneages_length, ended]=extraction(trajectories,DOX,min_length_traj)

%OUTPUT
%cycles_length : vector ça contains the durations of the cycles
%indeces : vector of the same length of cycles_length, containing repeted
%integers, for example k times the number i for the i-th lineage
%ligneages_length : vector contening the length of each lineage in
%generations
%ended : vector with length the number of lineages
%ended(j) = 1 if the j-th lineage is followed by a death, 0 otherwise

%INPUT
%trajectories : dataset (a struct)
%DOX : 1 for telomerase negative lineages only, 0 otherwise
%min_length_traj : minimal numebr of cycles in a lineages to include this
%lineages in the studied dataset

num_exp=numel(trajectories);
ended = zeros(1, numel(trajectories));
cycles_length=[];
indeces=[];
ligneages_length=[];

if DOX==0
    for i=1:num_exp
        if numel(trajectories(i).track)==0 %if empty trajectory
        continue
        end
        if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
            disp(['Film ' num2str(i) ' trop court']);
        else
            cycles_length=[cycles_length;trajectories(i).track(:,2)-trajectories(i).track(:,1)];
            indeces=[indeces;i*ones(size(trajectories(i).track(:,1)))];
            ligneages_length=[ligneages_length;length(trajectories(i).track(:,1))];
            ended(i) = trajectories(i).endlineage;
        end
    end
else % Only cycles after DOX addition
    for i=1:num_exp
        if numel(trajectories(i).track)==0 %if empty trajectory
        continue
        end
        if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
            disp(['Film ' num2str(i) ' trop court']);
        else
            t=(trajectories(i).track(:,2)-trajectories(i).track(:,1)).*(trajectories(i).track(:,1)>trajectories(i).DOXaddition); % Only cycles starting strictly after DOX addition
            cycles_length=[cycles_length;t(t~=0)];
            indeces=[indeces;i*ones(size(t(t~=0)))];
            ligneages_length=[ligneages_length;length(t(t~=0))];
            ended(i) = trajectories(i).endlineage;
        end
    end
end
