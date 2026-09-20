function dM = myGradmc(node,k,hK,c)
x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);
if k == 1
   dM = [[0+0*x, 0+0*x];[1+0*x, 0+0*x]./hK;[0+0*x, 1+0*x]./hK];
elseif k == 2
   dM =[[0+0*x, 0+0*x];[1+0*x, 0+0*x]./hK;[0+0*x, 1+0*x]./hK;[2*(x-xK),0+0*x]./hK^2; [(y-yK),(x-xK)]./hK^2; [0+0*x, 2*(y-yK)]./hK^2]; 
else
    disp('We are implementing case k ≠2');
end
