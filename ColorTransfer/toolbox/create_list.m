function Liste_pixels=create_list(im,Nx,Ny)
Liste_pixels=zeros(Nx*Ny,5);




im1=im(:,:,1);
im2=im(:,:,2);
im3=im(:,:,3);

[Y,X] = meshgrid(linspace(0,1,Ny), linspace(0,1,Nx));
Liste_pixels(:,4)=X(:)*255.;
Liste_pixels(:,5)=Y(:)*255.;

Liste_pixels(:,1)=im1(:);
Liste_pixels(:,2)=im2(:);
Liste_pixels(:,3)=im3(:);