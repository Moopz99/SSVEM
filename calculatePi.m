function Pi = calculatePi(elemID,centroid,diameter,node_new,elem_new)
k=2; 
indexK2 = elem_new{elemID};
hK = diameter(elemID);
x = node_new(indexK2 ,1);y = node_new(indexK2 ,2);
D = myM([x,y],k,hK,centroid(elemID,:));
Pi = (D'*D)\D';
Pi = blkdiag(Pi,Pi);
end
