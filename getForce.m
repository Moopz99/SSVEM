function ff = getForce(edgeMap,node_new,elem,press,direction)
k = 2;
sumNode = size(node_new,1);
ff = zeros(sumNode*2,1);
sumP = size(press,1);
nip = k+2;
[x,w] = gaussInt(nip);
for n = 1:sumP
    p = zeros(k+1,1);
    elemID = press(n,1);
    faceID = press(n,2);
    value = press(n,3);
    index = elem{elemID};
    Nv = length(index);
    v1 = 1:Nv; v2 = [2:Nv,1]; 
    elem1 = [v1(:), v2(:)];
    faceNode = elem1(faceID,:);
    faceNodeID = index(faceNode);
    edgeNodeCoor = node_new(faceNodeID,:);
    L = node_new(faceNodeID(1:2),:);
    L = (L(1,:)-L(2,:));
    Normal = [L(2),-L(1)]/norm(L);
    edgeKey = sprintf('%d_%d', min(faceNodeID), max(faceNodeID));
      if isKey(edgeMap, edgeKey)
          midPoint = edgeMap(edgeKey);
      else
          disp('Please check again');
          disp(edgeKey);
          pause;          
      end
      midCoord=node_new(midPoint,:);
      edgeNodeCoor =[edgeNodeCoor(1,:);midCoord;edgeNodeCoor(2,:)];
      faceNodeIDk2=[faceNodeID(1),midPoint,faceNodeID(2)];
      for m = 1:nip
          N = fun(x(m),0,0,1,k+1);
         dN = dfunc(x(m),0,0,1,k+1);
          J = sqrt((dN'*edgeNodeCoor(:,1))^2+(dN'*edgeNodeCoor(:,2))^2);
          p = p+w(m)*N*value*J;
      end
      if strcmp(direction,'normal')
          ff(faceNodeIDk2) = ff(faceNodeIDk2)+p*Normal(1);
          ff(faceNodeIDk2+sumNode) = ff(faceNodeIDk2+sumNode)+p*Normal(2);
      elseif strcmp(direction,'x')
          ff(faceNodeIDk2) = ff(faceNodeIDk2)+p;
      elseif strcmp(direction,'y')
          ff(faceNodeIDk2+sumNode) = ff(faceNodeIDk2+sumNode)+p;
      end
end
ff = sparse(ff);
end