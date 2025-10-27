function[cutoffselect]=GNM_sel_zjl(filename)
[~,resall,coorall,bexp]=read_pdb_protein(filename);
Pmax_b=0;
GNM_B=[];
E_B=[];
n=length(bexp);

for cutoff=5:18
    n1=n;
    netmat=zeros(n1); % Network matrix
    for i=1:n1
        for j=i+1:n1
            dis = sqrt((coorall(i,1)-coorall(j,1))^2+(coorall(i,2)-coorall(j,2))^2+(coorall(i,3)-coorall(j,3))^2);
            if dis<=cutoff
                netmat(i,j)=-1;
                netmat(j,i)=netmat(i,j);
            end
        end
    end
    for i=1:n1
        netmat(i,i)=-sum(netmat(i,:));
    end
    [V,D]=eig(netmat);
    d=zeros(n1,1);
    for i=1:n1
        d(i)=D(i,i);
    end
    flu=zeros(n1,1); % Fluctuation
    for i=1:n1
        for j=2:n1
            flu(i)=flu(i)+V(i,j)^2/d(j);
        end
    end

    N=n;
    B=regress(bexp,[ones(N,1) flu]);
    for i=1:n
        Bfactor_bsANM1(i)=B(1)+B(2)*flu(i);
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
    disp(PCC_b)
end

end