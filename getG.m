function G = getG(c,hK,node,indexK1,mat)
xE=c(1);
yE=c(2);
E = mat.E; nu =mat.nu;
ps=mat.ps;
 if ps==1         
   C = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; 
 elseif ps == 2   
   C = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))]; 
 else
     disp('Please check the flat condition');
 end
nBasis = 12; 
G = zeros(nBasis, nBasis); 
elemNodeCoords = node(indexK1, :);
nNodes = size(elemNodeCoords, 1);
triangles = cell(nNodes, 1);
for i = 1:nNodes
    j = mod(i,nNodes) + 1;
    triangles{i} = [xE, yE; elemNodeCoords(i, :); elemNodeCoords(j, :)];
end
for t = 1:nNodes
    triCoords = triangles{t};
    A = polyarea(triCoords(:,1), triCoords(:,2));
    if A < 0
        disp('Area is negative');
    end
    GaussPoints = [1/6, 1/6, 2/3 
                   2/3, 1/6, 1/6 
                   1/6, 2/3, 1/6];
        weights = [1/3, 1/3, 1/3]; 
    for gp = 1:3
        L = GaussPoints(gp, :);  
        x = L(1)*triCoords(1,1) + L(2)*triCoords(2,1) + L(3)*triCoords(3,1);
        y = L(1)*triCoords(1,2) + L(2)*triCoords(2,2) + L(3)*triCoords(3,2); 
        xi = (x - xE) / hK;
        eta = (y - yE) / hK;        
        strains = computeStrains(xi, eta, hK);
        for alpha = 1:nBasis
            for beta = 1:nBasis   
                term = strains(:, alpha)' * C * strains(:, beta);              
                G(alpha, beta) = G(alpha, beta) + weights(gp) * A * term;
            end
        end
    end
end
    
end
