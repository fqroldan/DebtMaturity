% Prueba

function [b1max,b2max,F,F1] = NC(tita0,alpha, delta,ph,nu,tau,gr,BU)

close all;
% By assumption tita0 > 1+delta, alpha >= 0.5, nu*(1-tau)+BL>0
BL=(nu*(1-tau)/10)*(-1);
B = linspace(BL,BU,gr);
b1 = zeros(1,gr);
b2 = zeros(1,gr);
b1n = zeros(1,gr);
b2n = zeros(1,gr);
be = zeros(gr,gr);
ben = zeros(gr,gr);
r = zeros(gr,gr);
r1 = zeros(gr,gr);



i=1;
while i<=size(b1,2);
     
    b1(1,i)=(nu*(1-tau)+B(i))^(1/2);
    
    i=i+1;   
end

i=1;
while i<=size(b2,2);
     
    b2(1,i)=(nu*(1-tau)+B(i))^(1/2);
    
    i=i+1;   
end

i=1;
while i<=size(b1n,2);
     
    b1n(1,i)=(nu*(1-tau)+B(i))^(-1/2);
    
    i=i+1;   
end

i=1;
while i<=size(b2n,2);
     
    b2n(1,i)=(nu*(1-tau)+B(i))^(-1/2);
    
    i=i+1;   
end
%Computes the expected value of b
i=1;
while i<=size(be,1);
     j=1;
    while j<=size(be,2);
    be(i,j)=(1+delta)^(1/2)*b1(i)*ph+(1-delta)^(1/2)*b1(i)*(1-ph)+((1+delta)*alpha+(1-delta)*(1-alpha))^(1/2)*b2(j)*ph+((1-delta)*alpha+(1+delta)*(1-alpha))^(1/2)*b2(j)*(1-ph);
    j=j+1;
    end;
    i=i+1;   
end

i=1;
while i<=size(be,1);
     j=1;
    while j<=size(be,2);
    ben(i,j)=(1+delta)^(1/2)*b1n(i)*ph+(1-delta)^(1/2)*b1n(i)*(1-ph)+((1+delta)*alpha+(1-delta)*(1-alpha))^(1/2)*b2n(j)*ph+((1-delta)*alpha+(1+delta)*(1-alpha))^(1/2)*b2n(j)*(1-ph);
    j=j+1;
    end;
    i=i+1;   
end



i=1;
while i<=size(r,2);
    j=1;
    while j<=size(r,1);
     
    r(i,j)=-tita0*(nu*(1-tau))*(3-(2*nu*(1-tau)*(ben(i,j))/(be(i,j))))^(-1)-(1/2)*((be(i,j))^(2));
    r1(i,j)=(nu*(1-tau))*(3-(2*nu*(1-tau)*(ben(i,j))/(be(i,j))))^(-1);
    
    j=j+1;
    end;
    i=i+1;   
end
F=r;
F1=r1;
[max_num,max_idx] = max(r(:));

[X Y]=ind2sub(size(r),max_idx);

b1max=B(1,X);
b2max=B(1,Y);

%figure(1)
%plot (B,r(:,1), B,r(:,100)), legend('b1 for min b2','b1 for max b2');

%figure(2)
%plot (B,r(1,:), B,r(100,:)), legend('b2 for min b1','b2 for max b1');