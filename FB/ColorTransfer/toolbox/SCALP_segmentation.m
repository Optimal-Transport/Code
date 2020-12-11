function [segmentation cluster]= SCALP_segmentation(im,N)
    
addpath('./SCALP/');




Nx=size(im,1);
Ny=size(im,2);
m = 0.02;  % Compactness parameter



disp('SCALP Super-pixel segmentation');

segments =SCALP_v1(im,N,m);
B = zeros(Nx,Ny);


%% Display SUPERPIXELS
% for i=1:max(segments(:))
%     j = (segments == i);
%     bb = bwboundaries(j);
%     if ~isempty(bb)
%         for k=1:length(bb{1})
%             B(bb{1}(k,1),bb{1}(k,2)) = 1;
%         end
%     end
%     
% end
% 
% figure,
% imagesc(double(im)/255.*repmat(~B,[1 1 3])); 

segments=segments+1;

nb_clusters=max(segments(:));% N*N;

Liste_pixels=zeros(Nx*Ny,5);


im1=im(:,:,1);
im2=im(:,:,2);
im3=im(:,:,3);

    
    [Y, X] = meshgrid(linspace(0,1,Ny), linspace(0,1,Nx));
Liste_pixels(:,4)=X(:)*255.;   %pour avoir la meme echelle que les couleurs
Liste_pixels(:,5)=Y(:)*255.;

Liste_pixels(:,1)=im1(:);
Liste_pixels(:,2)=im2(:);
Liste_pixels(:,3)=im3(:);



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
%weight_cluster=double(cpt_cluster)/Nx/Ny;




%calculer la moyenne  de chaque cluster


cluster=zeros(nb_clusters,5);

k_cluster=0;
for k=1:nb_clusters
    if length(Liste_cluster{k})>0  %enlever les clusters vides pour avoir une formule generale apres
        k_cluster=k_cluster+1;
        cluster(k_cluster,:)=mean(Liste_cluster{k});    
        I=segments==k;
        segments(I)=k_cluster;
    end
end

cluster=cluster(1:k_cluster,:);
segmentation=zeros(Nx*Ny,1);
segmentation(:)=segments(:);