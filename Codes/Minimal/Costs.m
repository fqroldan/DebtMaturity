% Prueba

function [CLC,CLI] = Costs(delta1,delta2,alpha,nu,tau,BM)

close all;
% By assumption tita0 > 1+delta, alpha >= 0.5, nu*(1-tau)+BL>0
deltagrid = linspace(delta1,delta2,length(BM));
b1 = zeros(1,length(BM));
b2 = zeros(1,length(BM));
C1H=zeros(1,length(BM));
C1L=zeros(1,length(BM));
CLC=zeros(1,length(BM));
CLI=zeros(1,length(BM));


i=1;
while i<=size(b1,2);
     
    b1(1,i)=(nu*(1-tau)+BM(1,i))^(1/2);
    
    i=i+1;   
end

i=1;
while i<=size(b2,2);
     
    b2(1,i)=(nu*(1-tau)+BM(2,i))^(1/2);
    
    i=i+1;   
end

%Computes the costs of Lack of Commitment and Lack of Insurance
i=1;
while i<=size(BM,2);
     
    
    C1H(i)=(1+deltagrid(i))^(-1/2)*b1(i)*((1+deltagrid(i))^(1/2)*b1(i)^(1/2)+((1+deltagrid(i))*alpha+(1-deltagrid(i))*(1-alpha))^(1/2)*b2(i));
    
    i=i+1;   
end

i=1;
while i<=size(BM,2);
     
    C1L(i)=(1-deltagrid(i))^(-1/2)*b1(i)*((1-deltagrid(i))^(1/2)*b1(i)+((1-deltagrid(i))*alpha+(1+deltagrid(i))*(1-alpha))^(1/2)*b2(i));
     
    i=i+1;   
end

i=1;
while i<=size(BM,2);
     
    CLI(i)=(C1H(i)/C1L(i))-((1-deltagrid(i))/(1+deltagrid(i)))^(1/2);
     
    i=i+1;   
end

i=1;
while i<=size(BM,2); %Cost of Lack of Commitment for the High state
     
    CLC(i)=((b1(i))^(2)/(b2(i))^(2))*(((1+deltagrid(i))*alpha+(1-deltagrid(i))*(1-alpha))/(1+deltagrid(i)))^(1/2)-(((1+deltagrid(i))*alpha+(1-deltagrid(i))*(1-alpha))/(1+deltagrid(i)))^(1/2);
     
    i=i+1;   
end

figure(1)
plot (deltagrid,CLC, deltagrid,CLI), legend('Cost of Lack of Commitment','Cost of Lack of insurance');

figure(2)
plot (deltagrid,BM(1,:), deltagrid,BM(2,:)), legend('b1','b2');




%figure(1)
%plot (B,r(:,1), B,r(:,100)), legend('b1 for min b2','b1 for max b2');

%figure(2)
%plot (B,r(1,:), B,r(100,:)), legend('b2 for min b1','b2 for max b1');