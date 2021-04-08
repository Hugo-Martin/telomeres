function [dureesdiv,indice,longueur_lignees, lignee_terminee]=extraction(trajectories,DOX,min_length_traj)

%OUTPUT
%dureesdiv : vecteur contenant les durées des divisions pour le jeu de données trajectories
%indice : vecteur de même dimension que dureesdiv, contenant une répétition
%de k fois le nombre i pour la i-ème lignée, qui est de longueur k
%longueur_lignee : vecteur contenant les longueur de chaque lignée (en nombre de
%divisions)
%lignee_terminee : vecteur de la taille du nombre de lignées,
%lignee_terminee(j) = 1 si la j-ème lignée se termine par une mort, 0 sinon

%INPUT
%trajectories : jeu de données
%ended : 1 si on ne considère que les lignées finies, 0 sinon
%DOX : 1 si on ne veut que les télomérase négatives, 0 sinon
%min_length_traj : nombre minimum de divisions dans une lignées pour ne pas
%la considérer trop courte

num_exp=numel(trajectories);
lignee_terminee = zeros(1, numel(trajectories));
dureesdiv=[];
indice=[];
longueur_lignees=[];

if DOX==0 %si on prend tous les cycles
    for i=1:num_exp
        if numel(trajectories(i).track)==0 %if empty trajectory
        continue
        end
        if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
            disp(['Film ' num2str(i) ' trop court']);
        else
            dureesdiv=[dureesdiv;trajectories(i).track(:,2)-trajectories(i).track(:,1)];
            indice=[indice;i*ones(size(trajectories(i).track(:,1)))];
            longueur_lignees=[longueur_lignees;length(trajectories(i).track(:,1))];
            lignee_terminee(i) = trajectories(i).endlineage;
        end
    end
else % si on ne prend que les cycles après l'ajout de DOX
    for i=1:num_exp
        if numel(trajectories(i).track)==0 %if empty trajectory
        continue
        end
        if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
            disp(['Film ' num2str(i) ' trop court']);
        else
            t=(trajectories(i).track(:,2)-trajectories(i).track(:,1)).*(trajectories(i).track(:,1)>trajectories(i).DOXaddition); % on ne prend que les cycles commençant strictement après l'ajout de DOX
            dureesdiv=[dureesdiv;t(t~=0)];
            indice=[indice;i*ones(size(t(t~=0)))];
            longueur_lignees=[longueur_lignees;length(t(t~=0))];
            lignee_terminee(i) = trajectories(i).endlineage;
        end
    end
end

% if (ended==0) && (DOX==0)
%     for i=1:num_exp
%         if numel(trajectories(i).track)==0 %if empty trajectory
%         continue
%         end
%         if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
%             disp(['Film ' num2str(i) ' trop court']);
%         else
%             dureesdiv=[dureesdiv;trajectories(i).track(:,2)-trajectories(i).track(:,1)];
%             indice=[indice;i*ones(size(trajectories(i).track(:,1)))];
%             longueur_lignees=[longueur_lignees;length(trajectories(i).track(:,1))];
%         end
%     end
% elseif (ended==1) && (DOX==0)
%     for i=1:num_exp
%         if numel(trajectories(i).track)==0 %if empty trajectory
%         continue
%         end
%         if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
%             disp(['Film ' num2str(i) ' trop court']);
%         else
%             if trajectories(i).endlineage==1
%             dureesdiv=[dureesdiv;trajectories(i).track(:,2)-trajectories(i).track(:,1)];
%             indice=[indice;i*ones(size(trajectories(i).track(:,1)))];
%             longueur_lignees=[longueur_lignees;length(trajectories(i).track(:,1))];
%             end
%         end
%     end
% elseif (ended==0) && (DOX==1)
%     for i=1:num_exp
%         if numel(trajectories(i).track)==0 %if empty trajectory
%         continue
%         end
%         if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
%             disp(['Film ' num2str(i) ' trop court']);
%         else
%             t=(trajectories(i).track(:,2)-trajectories(i).track(:,1)).*(trajectories(i).track(:,2)>=trajectories(i).DOXaddition);
%             dureesdiv=[dureesdiv;t(t~=0)];
%             indice=[indice;i*ones(size(t(t~=0)))];
%             longueur_lignees=[longueur_lignees;length(t(t~=0))];
%         end
%     end
% elseif (ended==1) && (DOX==1)
%     for i=1:num_exp
%         if numel(trajectories(i).track)==0 %if empty trajectory
%         continue
%         end
%         if length(trajectories(i).track(:,1))<min_length_traj %number of division cycles in the trajectory/lineage, 3 by default
%             disp(['Film ' num2str(i) ' trop court']);
%         else
%             if trajectories(i).endlineage==1
%                 t=(trajectories(i).track(:,2)-trajectories(i).track(:,1)).*(trajectories(i).track(:,2)>=trajectories(i).DOXaddition);
%                 dureesdiv=[dureesdiv;t(t~=0)];
%                 indice=[indice;i*ones(size(t(t~=0)))];
%                 longueur_lignees=[longueur_lignees;length(t(t~=0))];
%             end
%         end
%     end
% end