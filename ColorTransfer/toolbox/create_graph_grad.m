function K0=create_graph_grad(W,nb_neighbors)

%W is the n x n adjacency matrix (assumed symetric, but it should work if
%it is not)

nb_points=size(W,1);

%W=W+eye(nb_points)*1000; %to not select the correlation of the clusters with themselves


%create the NL sparse operator K0
K=zeros(nb_points,nb_neighbors);
X=zeros(nb_points,nb_neighbors);

for i=1:nb_points,
    A=W(i,:);
   
    %A=A+0.1;   %pour pas que ca explose, peut etre mieux de mettre un exp(-||.||)
    %A=double(ones(1,nb_points))./A;
    [val, index]=sort(A);   %ne garder que les n meilleurs voisins
    A(index(1:nb_points-nb_neighbors)) = 0;
 
    [I,J,V]=find(A);  
    K(i,1:length(J))=V;
    X(i,1:length(J))=J;
end

N = nb_points;
I =[1:N 1:N];
for k=1:nb_neighbors
    J = [X(:,k)' 1:N];
    % Have the value W) on the diagonal
    Identity = J(1:N)==1:N;
    V = [K(:,k)' -K(:,k)'];
    V(Identity) = 2*V(Identity);
    M = J>0;
    K0{k}=sparse(I(M),J(M),V(M),N,N);
end
