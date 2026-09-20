function [stress,mises,smax,smin,strain] = calculateStress(node,elem,node_new,elem_new,uh,mat)
sumElem = size(elem,1); 
sumNode = size(node_new,1);
ux = uh(1:sumNode);
uy = uh(sumNode+1:end);
stress = zeros(sumNode,3);
strain = zeros(sumNode,3);
nodeUsed = zeros(sumNode,1);
centroid = zeros(sumElem,2); diameter = zeros(sumElem,1);
s = 1;
for iel = 1:sumElem
    indexK1 = elem{iel};
    verts = node(indexK1, :); verts1 = circshift(verts,-1);
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    centroid(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
    diameter(s) = max(pdist(verts));    
    s = s+1; 
end
 E = mat.E; nu =mat.nu; ps=mat.ps;
 if  ps==1     
 C = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
 elseif ps==2  
 C = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))]; 
 else
     disp('calculateStress error');
 end
for n = 1:sumElem
    A1 = [1,0;0,0;0,1];
    A2 = [0,0;0,1;1,0];
    A = [A1,A2];
    Pi = calculatePi(n,centroid,diameter,node_new,elem_new);  
    indexK2 =elem_new{n};
    NvK2 = length(indexK2);
    x=node_new(indexK2,1); y=node_new(indexK2,2);
    xy=[x,y];
    for i=1:NvK2     
    dM=myGradmc(xy(i,:),2,diameter(n),centroid(n,:));
    dM1 = blkdiag(dM',dM');
    ee = A*dM1*Pi*[ux(indexK2);uy(indexK2)];
    eS = C*ee;
    strain(indexK2(i),:) = strain(indexK2(i),:)+ee';
    stress(indexK2(i),:) = stress(indexK2(i),:)+eS';
    nodeUsed(indexK2(i)) = nodeUsed(indexK2(i))+1;
    end
end
strain = strain./nodeUsed;
stress = stress./nodeUsed;
smax = (stress(:,1)+stress(:,2))./2+sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);
smin = (stress(:,1)+stress(:,2))./2-sqrt(((stress(:,1)-stress(:,2))./2).^2+stress(:,3).^2);
mises = sqrt((smax.^2+smin.^2+(smax-smin).^2)./2);
end
