function D = myD_matrix(node,k,hK,c)
x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);
xi = (x-xK)./hK; 
eta = (y-yK)./hK;
if k == 2 
D = [[1.0+0*x;0*x],[0*x;1.0+0*x],[xi;0*x],[0*x;xi],[eta;0*x],[0*x;eta],...
    [xi.*xi;0*x],[0*x;xi.*xi],[xi.*eta;0*x],[0*x;xi.*eta],[eta.*eta;0*x],[0*x;eta.*eta]];
else
disp('We are implementing case k ≠2');
end

