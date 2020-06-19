% Prueba

function [BM] = mat(delta1,delta2,gr,tita0,alpha,ph,nu,tau,BU)

close all;

deltagrid = linspace(delta1,delta2,gr);
r = zeros(2,max(size(deltagrid)));

i=1;
while i<=size(r,2);
        [b1max,b2max] = NC1(tita0,alpha, deltagrid(i),ph,nu,tau,gr,BU);
    r(1,i)=b1max;
    r(2,i)=b2max;
     
    i=i+1;   
end
BM=r;
figure(1)
plot (deltagrid,r(1,:), deltagrid,r(2,:)), legend('b1','b2');
