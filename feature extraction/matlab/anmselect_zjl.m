function [cutoffselect]=anmselect_zjl(filename)
Pmax_b=0;
cutoffselect=0;
GNM_B=[];
E_B=[];
[~,resall,coorall,bexp]=read_pdb_protein(filename);
n=length(resall);
l = [];
s=[];
hessian=zeros;

for cutoff=7:20 % Test cutoff range 7:20, find cutoff with maximum PCC
    hessian=zeros(3*n);
    for i=1:n
        for j=i+1:n
            dis=sqrt(sum((coorall(i,:)-coorall(j,:)).^2));
            if dis == 0
                l = [l;i];
                coorall(i,:)
                coorall(j,:)
            end
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
    fluslown1=zeros(n,1);
    for i=1:n
        for k=7:3*n
            fluslown1(i)=fluslown1(i)+V((3*i-2),k)*V((3*i-2),k)/D(k,k)+V((3*i-1),k)*V((3*i-1),k)/D(k,k)+V((3*i),k)*V((3*i),k)/D(k,k);
        end
    end
    
    N=n;
    B=regress(bexp,[ones(N,1) fluslown1]);
    for i=1:n
        Bfactor_bsANM1(i)=B(1)+B(2)*fluslown1(i);
    end
    [S3]=corrcoef(bexp,Bfactor_bsANM1);
    PCC_b(cutoff)=S3(1,2);
    E_B=bexp;
    dfactor(:,2)=Bfactor_bsANM1;
    dfactor(:,1)=E_B;
    [S,P]=corrcoef(dfactor); % S is correlation coefficient n2x2 matrix
    if S(1,2)>=Pmax_b
        Pmax_b=S(1,2);
        cutoffselect=cutoff;
    end
end

end