function [Im_final]=synthesis(Liste_pixels,Moyenne_clusterX,Variance_clusterX,Moyenne_clusterY,Flow,Nx,Ny,nb_cluster)


Im_final=zeros(Nx*Ny,3);
poids_final=zeros(Nx*Ny,1);
Moyenne_clusterXY=zeros(nb_cluster,5);
A=zeros(Nx*Ny,5);
for k=1:nb_cluster
    
    %creation des couleurs interpolees
    Moyenne_clusterXY(k,:)=(Flow(k,:)*Moyenne_clusterY)/(sum(Flow(k,:)));

    
    %calcul de la vraissemblance des pixels aux clusters d'origine
    for l=1:5
        A(:,l)=Liste_pixels(:,l)-Moyenne_clusterX(k,l);
    end
    
    % debug Julien : ajout /255 -> mettre 255^2 ?
    A = A / 255; 
    
    V = Variance_clusterX{k}; % Attention : si variance trop grande, les poids deviennent nuls ... problème à régler
    poids_tmp=exp(-sum((A*V).*A,2));

    poids_final=poids_final+poids_tmp;
    
    for l=1:3
        Im_final(:,l)=Im_final(:,l)+poids_tmp*Moyenne_clusterXY(k,l);
    end
    
end

for l=1:3
    Im_final(:,l)=Im_final(:,l)./poids_final(:);
end


Im_final=Im_final';
