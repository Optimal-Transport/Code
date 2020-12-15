function [Moyenne_cluster Variance_cluster]=compute_stat_clusters(Liste_pixels,segments,nb_clusters,Nx,Ny,sigma_color,sigma_spatial)

%recuperer les pixels de chaque cluster
for k=1:nb_clusters
    Liste_cluster{k}=[];
end


cpt_cluster=zeros(nb_clusters,1);

for i=1:Nx*Ny,
    k=segments(i);
    cpt_cluster(k)=cpt_cluster(k)+1;
    Liste_cluster{k}(cpt_cluster(k),:)=Liste_pixels(i,:);
    
end

%poids du cluster



%calculer la moyenne et la variance de chaque cluster

Filtre_variance=ones(5,5);
%

Filtre_variance=[1. 1. 1. 0. 0.;1. 1. 1. 0. 0.;  1. 1. 1. 0 0;  0 0 0 1. 1.;  0 0 0 1. 1.];
%pour decoreller spatial et colorimetrique
%Filtre_variance=eye(5);
scal=ones(5,5);%*2.-eye(5);   
%on ne se servira que de la matrice tri sup apres pour calculer les vraissemblances, donc on double les valeurs de la covariance en dehors de la diagonale


vect=[sigma_color sigma_color sigma_color sigma_spatial sigma_spatial ];
Filtre_variance=(vect'*vect).*Filtre_variance;

Epsilon=0.0001*eye(5);
Moyenne_cluster=ones(nb_clusters,5);


for k=1:nb_clusters
%    Liste_cluster{k}
    if length(Liste_cluster{k})==0  %mettre infini si cluster vide, impossible normalement
        Moyenne_cluster(k,:)=rand(1,5)*1000000.;
        disp('Problem');
        Variance_cluster{k}=ones(5,5)*1000.;
    else
        Moyenne_cluster(k,:)=mean(Liste_cluster{k});
 

  
         Variance_cluster{k}=(inv(cov(Liste_cluster{k}/255.).*Filtre_variance+Epsilon).*scal);
         
%         Variance_cluster{k}=zeros(5,5);
%         Variance_cluster{k}(1:3,1:3)=(inv(cov(Liste_cluster{k}(:,1:3)/255.).*Filtre_variance(1:3,1:3)+Epsilon(1:3,1:3)).*scal(1:3,1:3));
%         Variance_cluster{k}(4:5,4:5)=(inv(cov(Liste_cluster{k}(:,4:5)/255.).*Filtre_variance(4:5,4:5)+Epsilon(4:5,4:5)).*scal(4:5,4:5));
    end
end