function strains = computeStrains(xi, eta, hK)
strains = zeros(3, 12); 
strains(:, 1) = [0; 0; 0];        
strains(:, 2) = [0; 0; 0];        
strains(:, 3) = [1/hK; 0; 0];      
strains(:, 4) = [0; 0; 1/hK];    
strains(:, 5) = [0; 0; 1/hK];     
strains(:, 6) = [0; 1/hK; 0];     
strains(:, 7) = [2*xi/hK; 0; 0]; 
strains(:, 8) = [0; 0; 2*xi/hK];  
strains(:, 9) = [eta/hK; 0; xi/hK];  
strains(:, 10) = [0; xi/hK; eta/hK]; 
strains(:, 11) = [0; 0; 2*eta/hK];   
strains(:, 12) = [0; 2*eta/hK; 0];   
end