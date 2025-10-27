function[col,b]=ANM_calPRS_zjl(filename)
Pmax_b=0;
GNM_B=[];
E_B=[];
[~,resall,coorall,bexp]=read_pdb_protein(filename);
n=length(resall);
s=[];
hessian=zeros;
c = anmselect_zjl(filename); 

%% Build network
for cutoff=c
    hessian=zeros(3*n);
    for i=1:n
        for j=(i+1):n
            dis=sqrt(sum((coorall(i,:)-coorall(j,:)).^2));
            if (dis<=cutoff)
                hessian(3*i-2,3*j-2)=-(coorall(j,1)-coorall(i,1))*(coorall(j,1)-coorall(i,1))/dis^2;
                hessian(3*i-2,3*j-1)=-(coorall(j,1)-coorall(i,1))*(coorall(j,2)-coorall(i,2))/dis^2;
                hessian(3*i-2,3*j)=-(coorall(j,1)-coorall(i,1))*(coorall(j,3)-coorall(i,3))/dis^2;
                hessian(3*i-1,3*j-2)=hessian(3*i-2,3*j-1);
                hessian(3*i-1,3*j-1)=-(coorall(j,2)-coorall(i,2))*(coorall(j,2)-coorall(i,2))/dis^2;
                hessian(3*i-1,3*j)=-(coorall(j,2)-coorall(i,2))*(coorall(j,3)-coorall(i,3))/dis^2;
                hessian(3*i,3*j-2)=hessian(3*i-2,3*j);
                hessian(3*i,3*j-1)=hessian(3*i-1,3*j);
                hessian(3*i,3*j)=-(coorall(j,3)-coorall(i,3))*(coorall(j,3)-coorall(i,3))/dis^2;
            end
            hessian(3*j-2,3*i-2)=hessian(3*i-2,3*j-2);
            hessian(3*j-1,3*i-2)=hessian(3*i-2,3*j-1);
            hessian(3*j,3*i-2)=hessian(3*i-2,3*j);
            hessian(3*j-2,3*i-1)=hessian(3*i-1,3*j-2);
            hessian(3*j-1,3*i-1)=hessian(3*i-1,3*j-1);
            hessian(3*j,3*i-1)=hessian(3*i-1,3*j);
            hessian(3*j-2,3*i)=hessian(3*i,3*j-2);
            hessian(3*j-1,3*i)=hessian(3*i,3*j-1);
            hessian(3*j,3*i)=hessian(3*i,3*j);
        end
    end
    for i=1:n
        for j=1:n
            dis=sqrt(sum((coorall(i,:)-coorall(j,:)).^2));
            if ((j~=i) & (dis<=cutoff))
                hessian(3*i-2,3*i-2)=hessian(3*i-2,3*i-2)-hessian(3*i-2,3*j-2);
                hessian(3*i-2,3*i-1)=hessian(3*i-2,3*i-1)-hessian(3*i-2,3*j-1);
                hessian(3*i-2,3*i)=hessian(3*i-2,3*i)-hessian(3*i-2,3*j);
                hessian(3*i-1,3*i-2)=hessian(3*i-2,3*i-1);
                hessian(3*i-1,3*i-1)=hessian(3*i-1,3*i-1)-hessian(3*i-1,3*j-1);
                hessian(3*i-1,3*i)=hessian(3*i-1,3*i)-hessian(3*i-1,3*j);
                hessian(3*i,3*i-2)=hessian(3*i-2,3*i);
                hessian(3*i,3*i-1)=hessian(3*i-1,3*i);
                hessian(3*i,3*i)=hessian(3*i,3*i)-hessian(3*i,3*j);
            end
        end
    end
    [V,D]=eig(hessian); 
    bfactors=zeros(n,1);
end % Network building complete

U=V';
hessian_inverse=zeros(3*n);
for k=7:3*n
    hessian_inverse=hessian_inverse+V(:,k)*U(k,:)/D(k,k);
end

%% Calculate PRS matrices
R1=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=1;
    F(3*j-1)=1;
    F(3*j)=1;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R1=R1+dR;
end

R2=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=1;
    F(3*j-1)=1;
    F(3*j)=0;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R2=R2+dR;
end

R3=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=1;
    F(3*j-1)=0;
    F(3*j)=1;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R3=R3+dR;
end

R4=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=0;
    F(3*j-1)=1;
    F(3*j)=1;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R4=R4+dR;
end

R5=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=1;
    F(3*j-1)=0;
    F(3*j)=0;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R5=R5+dR;
end

R6=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=0;
    F(3*j-1)=1;
    F(3*j)=0;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R6=R6+dR;
end

R7=zeros(n);
for j=1:n
    F=zeros(3*n,1);
    dR=zeros(n);
    F(3*j-2)=0;
    F(3*j-1)=0;
    F(3*j)=1;
    ddR=hessian_inverse*F;
    W=zeros(n,1);
    for i=1:n
        W(i)=sqrt(ddR(3*i-2)*ddR(3*i-2)+ddR(3*i-1)*ddR(3*i-1)+ddR(3*i)*ddR(3*i));
    end
    dR(:,j)=W;
    R7=R7+dR;
end

R=(R1+R2+R3+R4+R5+R6+R7)/7;
W=diag(R);
M=1./W;
N=diag(M);
P=N*R;
figure(1)
map=jet(n);
colormap(map);
pcolor(P);
res=[1,n,1,n];
axis(res)
shading interp;
view(0,90);
clim([0 1]);
colorbar;
title('PRS');
a=mean(P,1); % Column average
PRS2=a';
b=mean(P,2); % Row average
x= 1:1:n;
figure(2)
plot(x,a);
title('PRS Column');
figure(3)
plot(x,b);
title('PRS Line');
col = a';

end