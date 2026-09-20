function GK = globalK(node,elem,node_new,elem_new,mat)
sumElem = size(elem,1); 
sumNode = size(node_new,1);
centroid = zeros(sumElem,2); diameter = zeros(sumElem,1);
s = 1;
for iel = 1:sumElem
    index = elem{iel};
    verts = node(index, :); verts1 = circshift(verts,-1);
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    centroid(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
    diameter(s) = max(pdist(verts));    
    s = s+1;
end
elemLen = cellfun('length',elem_new);
nnz = sum((2*elemLen).^2); 
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
ia = 0;
for n = 1:sumElem
    AK = elemK(n,centroid,diameter,node_new,elem_new,node,elem,mat);
    AB = reshape(AK',1,[]);
    index = elem_new{n};
    indexDof = [index, index+sumNode]; 
    Ndof = length(indexDof);
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
end
GK = sparse(ii,jj,ss,sumNode*2,sumNode*2);
