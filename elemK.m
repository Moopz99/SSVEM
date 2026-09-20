function K = elemK(elemID,centroid,diameter,node_new,elem_new,node,elem,mat)
k=2; 
a = size(elem,1);
b = size(node_new,1);
indexK2 = elem_new{elemID};
hK = diameter(elemID);
indexK1 = elem{elemID};
x = node_new(indexK2 ,1);y = node_new(indexK2 ,2);
D = myD_matrix([x,y],k,hK,centroid(elemID,:)); 
Pi = (D'*D)\D';
G = getG(centroid(elemID,:),hK,node,indexK1,mat);
KEc = Pi'*G*Pi;
DPi=D*Pi;
I = eye(size(DPi));
alpha0=1.0;
Nen=size(KEc,1);
alpha=0.0;
for i=1:Nen
    alpha=alpha+KEc(i,i);
end
alpha=(alpha0/(1.0*Nen))*(a/b)*alpha;
KEs = alpha*((I-DPi)' *(I - DPi));
K=KEc+KEs;
end

