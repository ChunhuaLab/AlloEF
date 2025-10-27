function [bet,zzz]=comnetbet(filename)
[~,resall,coorall,bexp]=read_pdb_protein(filename);
n=length(resall);
netmat=zeros(n);

for i=1:n
    for j=i:n
        if i==j
            continue;
        else
            x1=coorall(i,:);
            x2 = coorall(j,:);
            dis=sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2+(x1(3)-x2(3))^2);
            cutoff =7; % Cutoff radius
            if dis<cutoff
                netmat(i,j)=1;
                netmat(j,i)=1;
            end
        end
    end
end

for i=1:n
    netmat(i,i)=sum(netmat(i,:));
end
avenetmat=diag(netmat);

%% Calculate network degree
zzz=zeros(n,1);
for i=1:n
    for j=1:n
        if netmat(i,j)~=0
            zzz(i,1)=zzz(i,1)+netmat(i,j);
        end
    end
end
% figure(1)
% plot(zzz)
% title('Node Degree');
% xlabel('Residue');
% ylabel('Node Degree');
% hold on

%% Calculate shortest path
C=zeros(n,n);
Dis=netmat;
C=netmat;
for i =1:n
    for j=1:n
        if Dis(i,j)==0 && i~=j
            Dis(i,j)=inf;
        end
    end
end

for k=1:n % Floyd algorithm to find shortest path length between any two points
    for i=1:n
        for j=1:n % Can calculate only half since it's symmetric, doesn't affect result size, but remember to add D(j,i)=D(i,j)
            if Dis(i,j)>Dis(i,k)+Dis(k,j)
                Dis(i,j)=Dis(i,k)+Dis(k,j); % Update distance between i and j
                C(i,j)=C(i,k)*C(k,j); % Update shortest path count
            elseif Dis(i,j)==Dis(i,k)+Dis(k,j)
                if k~=i&&k~=j % To avoid double counting, exclude endpoints here
                    C(i,j)=C(i,j)+C(i,k)*C(k,j); % Update path count with same shortest distance
                end
            end
        end
    end
end

figure(3);
m1=64;
map=jet(m1);
colormap(map);
pcolor(Dis);
axis square;
xlabel('Residue');
ylabel('Residue');
shading interp;
view(0,90);
colorbar('FontSize',15,'FontName','Times New Roman');
title('Shortest Path');

%% Calculate betweenness
B1=zeros(1,n);
for k=1:n
    for i=1:n
        if i~=k % Exclude endpoints
            for j=i+1:n % Since it's undirected, calculate only forward direction once, so only half; if calculating both directions, sum would be double, doesn't affect relative values
                if j~=k % Exclude endpoints
                    if Dis(i,j)==Dis(i,k)+Dis(k,j)&C(i,j)~=0 % Condition proves shortest path between i and j passes through k node
                        B1(k)=B1(k)+C(i,k)*C(k,j)/C(i,j);
                    end
                end
            end
        end
    end
end

figure(5)
plot(B1)
title('Betweenness');
xlabel('Residue');
ylabel('Betweenness');

%% Normalize
zbzh=zscore(B1);
bet = zbzh';
% figure(6)
% plot(zbzh)
% xlabel('Residue');
% ylabel('Betweenness');
% title('Normalized Betweenness')

end