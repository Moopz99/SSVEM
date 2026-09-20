function D = myM(node,k,hK,c)
x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);
xi = (x-xK)./hK; 
eta = (y-yK)./hK;
if k == 2 
D = [1.0+0*x,xi,eta,xi.*xi,xi.*eta,eta.*eta];
else
disp('We are implementing case k ≠2');
end