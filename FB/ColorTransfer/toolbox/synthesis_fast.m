function Im_final=synthesis_fast(Liste_pixels,Moyenne_clusterX,Variance_clusterX,Moyenne_clusterY,Flow,Nx,Ny,nb_cluster)

Im_final=zeros(Nx*Ny,3);
poids_final=zeros(Nx*Ny,1);
A=zeros(Nx*Ny,5);
MoyenneXY=diag(1./sum(Flow,2))*(Flow(:,:)*Moyenne_clusterY(:,1:3));%

%Vectorialisation

Pixels=repmat(Liste_pixels,[1 nb_cluster]);

Clusters=repmat(reshape(Moyenne_clusterX',1,nb_cluster*5),[Nx*Ny 1]);


moyenne=(Pixels-Clusters)/255.;

clear  Pixels Clusters

poids=moyenne*Variance_clusterX;

poids=poids.*moyenne;
% toc
% poids=Variance_clusterX*moyenne';
% poids=poids'.*moyenne;

clear moyenne;

poids=reshape(exp(-sum(reshape(poids,Nx*Ny,5,nb_cluster),2)),[Nx*Ny nb_cluster]);
poids2=repmat(sum(poids,2), [1 3]);
poids2(poids2==0)=1.;

Im_final=(poids*MoyenneXY./ poids2)';



% 
% 
% for k=1:nb_cluster
%     %creation des couleurs interpolees
%     %Moyenne_clusterXY=MoyenneXY(k,:)/sumXY(k);
%    
%  
%     %calcul de la vraissemblance des pixels aux clusters d'origine
%      poids_tmp = 0;
%       for l=1:5
%           A(:,l) = (Liste_pixels(:,l)-Moyenne_clusterX(k,l))/255.;
%           poids_tmp=poids_tmp+A(:,l).^2*Variance_clusterX{k}(l,l);
%       end
% %  
% 
%     %  A=(Liste_pixels-repmat(Moyenne_clusterX(k,:), [Nx*Ny 1]))/255.;
% 
%     %  poids_tmp=exp(-sum((A*Variance_clusterX{k}).*A,2));
%     
%     % hypothèse : Matrice de coVariance symmetrique
%     % Julien : il ne faudrait pas ajouter un facteur '2' du coup ?
%     %si mais c'est dans la variance, cf ligne 31 de compute_stat_clusters
%     for l=1:4
%         for l2=l+1:5
%             poids_tmp = poids_tmp +  A(:,l) .* A(:,l2) * Variance_clusterX{k}(l,l2);
%         end
%     end
%     
%     % exponential weight kernel (gaussian mixture model)
%     poids_tmp=exp(-poids_tmp);
%   
%     poids_final = poids_final + poids_tmp;
%    
%     Im_final = Im_final + poids_tmp * MoyenneXY(k,:);
%     
% end
% 
% Im_final = Im_final' ./ repmat(poids_final', [3 1]);
%toc
% max(Im_final(:)-Im_final0(:))
