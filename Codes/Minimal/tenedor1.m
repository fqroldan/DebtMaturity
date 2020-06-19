% Prueba

function tenedor1

close all;
%tenedor(alpha, delta,nu,tau,gr,BU,psi,ph)
ph=0.5;
psi=0.5;% FRAN ACÁ LO TENES QUE CAMBIAR A MANO
BU=0.5;
gr=100; % FRAN ACÁ LO TENES QUE CAMBIAR A MANO (pone 400 por ejemplo, es el tamaño de la grilla).
tau=0.1;
nu=1;
delta=0.3; % FRAN ACÁ LO TENES QUE CAMBIAR A MANO
alpha=0.7;
% By assumption tita0 > 1+delta, alpha >= 0.5, nu*(1-tau)+BL>0
BL=0;
B = linspace(BL,BU,gr);
B02 = linspace(BL,BU,gr);
B01 = linspace(BL,BU,gr);
beta=0.99;
tita0=2;

c2 = zeros(1,gr);
g2 = zeros(1,gr);
x1 = zeros(gr,gr);
c1 = cell(1,gr);
tita1H=1+delta;
tita1L=1-delta;
tita2H=alpha*(1+delta)+(1-alpha)*(1-delta);
tita2L=alpha*(1-delta)+(1-alpha)*(1+delta);
OF1H= cell(1,gr);
OF1L= cell(1,gr);
OF2= zeros(gr,gr);
EC= zeros(gr,gr);
C0= zeros(gr,gr);
%OF2= zeros(gr,gr);
B12Hmax= zeros(gr,gr);
B12Lmax= zeros(gr,gr);




i=1;
while i<=size(g2,2); % i is B12
     
    g2(1,i)=tau*nu-B(1,i);
    
    i=i+1;   
end

i=1;
while i<=size(c2,2); % i is B12
     
    c2(1,i)=nu-g2(1,i);
    
    i=i+1;   
end

i=1;
while i<=size(x1,1);% i is B12
     j=1;% j is B02
    while j<=size(x1,2);
    x1(i,j)=(c2(1,i)-nu*(1-tau)-B02(1,j))/c2(1,i);
    j=j+1;
    end;
    i=i+1;   
end

i=1;
while i<=gr;% i is B12
   j=1;% j is B02
   while j<=gr;
      h=1;% h is B01
    while h<=gr;
  
    c1{1,h}(i,j)=(B01(1,h)+nu*(1-tau))/(1+x1(i,j));
    
    h=h+1;
    end
    j=j+1;
  end;
    i=i+1;   
end


%Computes the expected value of b KEEP ON HERE
i=1;
while i<=gr;% i is B12
   j=1;% j is B02
   while j<=gr;
      h=1;% h is B01
    while h<=gr;

        OF1H{1,h}(i,j)=psi*(tita1H*(nu-c1{1,h}(i,j))+beta*tita2H*g2(1,i))+(1-psi)*(log(c1{1,h}(i,j))+beta*log(c2(1,i)));
    
        h=h+1;
    end
    j=j+1;
  end;
    i=i+1;   
end

i=1;
while i<=gr;% i is B12
   j=1;% j is B02
   while j<=gr;
      h=1;% h is B01
    while h<=gr;

        OF1L{1,h}(i,j)=psi*(tita1L*(nu-c1{1,h}(i,j))+beta*tita2L*g2(1,i))+(1-psi)*(log(c1{1,h}(i,j))+beta*log(c2(1,i)));
    
        h=h+1;
    end
    j=j+1;
  end;
    i=i+1;   
end



for j=1:gr % j is B02
    for h=1:gr % h is B01
[max_num,max_idx] = max(OF1H{1,h}(:,j));
B12Hmax(j,h)=max_idx;
    end
end

for j=1:gr % j is B02
    for h=1:gr % h is B01
[max_num,max_idx] = max(OF1L{1,h}(:,j));
B12Lmax(j,h)=max_idx;
    end
end
b12H=B12Hmax;
b12L=B12Lmax;

for j=1:gr % j is B02
    for h=1:gr % h is B01
EC(j,h)=ph*((c1{1,h}(b12H(j,h),j)-nu*(1-tau))/c1{1,h}(b12H(j,h),j))+(1-ph)*((c1{1,h}(b12L(j,h),j)-nu*(1-tau))/c1{1,h}(b12L(j,h),j))+ph*((c2(1,b12H(j,h))-nu*(1-tau))/c2(1,b12H(j,h)))+(1-ph)*((c2(1,b12L(j,h))-nu*(1-tau))/c2(1,b12L(j,h)));
    end
end

for j=1:gr % j is B02
    for h=1:gr % h is B01
C0(j,h)=(nu*(1-tau)/(EC(j,h)+1));
    end
end

for j=1:gr % j is B02
    for h=1:gr % h is B01
OF2(j,h)=(1-psi)*(log(C0(j,h))+beta*(ph*log(c1{1,h}(b12H(j,h),j))+(1-ph)*log(c1{1,h}(b12L(j,h),j)))+beta*beta*(ph*log(c2(1,b12H(j,h)))+(1-ph)*log(c2(1,b12H(j,h)))))+psi*(tita0*(nu-C0(j,h))+beta*(ph*tita1H*(nu-c1{1,h}(b12H(j,h),j))+(1-ph)*tita1L*(nu-c1{1,h}(b12L(j,h),j)))+beta*beta*(ph*tita2H*(nu-c2(1,b12H(j,h)))+(1-ph)*tita2L*(nu-c2(1,b12H(j,h)))));
    end
end

[max_num1,max_idx1] = max(OF2(:));
[B02OPT B01OPT]=ind2sub(size(OF2),max_idx1);
bhat02=B(1,B02OPT);
bhat01=B(1,B01OPT);

save('tenedor1.mat','bhat02','bhat01')
end