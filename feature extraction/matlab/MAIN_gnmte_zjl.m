clear;
clc;

%% 
warningFileName = '.\protein_warnings.txt';

stateFileName = 'processing_state.mat';

% Load processing state if exists
if exist(stateFileName, 'file')
    load(stateFileName, 'startIdx'); % Resume from previous state
    fprintf('Resuming processing from file index %d\n', startIdx);
else
    startIdx = 1; % Start from beginning
end

fid = fopen(warningFileName, 'w');
folderPath = '.\1WQW_A.pdb'; % Folder path
files = dir(fullfile(folderPath, '*.pdb')); % Get all pdb files in folder

for fileIdx = 1:length(files)
    filename = fullfile(folderPath, files(fileIdx).name);
    [residall, resall, coorall, bexp,chainall] = read_pdb_protein(filename);
    [n, ~] = size(coorall);
    
    %     if n > 2200
    %     fprintf(fid, 'File %s: Atom count greater than 3000, skipping this file\n', files(fileIdx).name);
    %     continue; % Skip processing
    % end

    if n > 2000
        fprintf(fid, 'File %s: Atom count greater than 2000\n', files(fileIdx).name);
    end

    cutoff = GNM_sel_zjl(filename); % Use custom defined cut off: optimized cutoff

    % Build network
    Kirchhoff = zeros(n); % Connection matrix
    for i = 1:n
        for j = i+1:n
            dis = sqrt(sum((coorall(i,:) - coorall(j,:)).^2));
            if dis <= cutoff
                Kirchhoff(i,j) = -1;
                Kirchhoff(j,i) = Kirchhoff(i,j);
            end
        end
    end

    for i = 1:n
        Kirchhoff(i,i) = -sum(Kirchhoff(i,:));
    end

    [V,D] = eig(Kirchhoff);

    % You can add result processing and saving here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate msf
n1=n;
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
PCC_b=S3(1,2);

if PCC_b < 0.5
    fprintf(fid, 'File %s: PCC_b less than 0.5\n', files(fileIdx).name);
end

%% calculate PRS
[col_n,line_n]=ANM_calPRS_zjl(filename);

%% transfer entropy
[norm_ave5, netTE5, norm_netTE5] = Nte(5, n, V, D);
[norm_ave15, netTE15, norm_netTE15] = Nte(15, n, V, D);

saveFolderPath = 'C:\Users\13017\Desktop\allo protein prediction\feature\independent_test_feature\NTE-sites\result';

% Get the file name without extension
[~, baseFilename, ~] = fileparts(filename);

% Create filenames for the Excel files
netTEFileName = fullfile(saveFolderPath, strcat(baseFilename, '_netTE.xlsx'));
normNetTEFileName = fullfile(saveFolderPath, strcat(baseFilename, '_norm_netTE.xlsx'));

% Adding headers and row names to netTE matrix, including chain name
resNames = strcat(resall, ' (', string(residall), ')', ' - ', chainall); % Create the name-sequence-chain labels

% Find unique indices for resNames and remove duplicates
[uniqueResNames, uniqueIndices] = unique(resNames, 'stable');

% Use only unique indices for creating the tables
netTE_table = array2table(netTE15(uniqueIndices, uniqueIndices), 'VariableNames', uniqueResNames, 'RowNames', uniqueResNames);
norm_netTE_table = array2table(norm_netTE15(uniqueIndices, uniqueIndices), 'VariableNames', uniqueResNames, 'RowNames', uniqueResNames);

% Save matrices to Excel files
writetable(netTE_table, netTEFileName, 'WriteRowNames', true);
writetable(norm_netTE_table, normNetTEFileName, 'WriteRowNames', true);

%% calculate betweenness
[norm_bet,degree]=comnetbet(filename);

% Output results
[~, filename, ext] = fileparts(filename);
filenameWithExt = strcat(filename, ext);
resall_cell = cellstr(resall); % Convert character matrix to cell array
outputTable = table(residall, resall_cell, flu,norm_ave5,norm_ave15, col_n,line_n,norm_bet,degree,chainall,...
    'VariableNames', {'ResidueIndex', 'ResidueInfo', 'msf','NTE5','NTE15','PRScol','PRSlin','norm_bet','degree','chain'});

filePath = sprintf('C:\\Users\\13017\\Desktop\\allo protein prediction\\feature\\independent_test_feature\\matlab\\result\\%s.csv',filenameWithExt);

% Save table to CSV file at specified location
writetable(outputTable, filePath);

%% Clear data
clear resall coorall bexp;
clear Kirchhoff V D B;
clear flu Bfactor_bsANM1 S3 col_n line_n norm_ave5 norm_ave15 norm_bet degree;

if fid == -1
    error('File could not be opened.');
end
end

fclose(fid);
if exist(stateFileName, 'file')
    delete(stateFileName);
end
disp('Processing completed!');

%% transfer formula calculation

function [norm_ave,netTE,norm_netTE]=Nte(tau,n,V,D)

% tau=15; Recommended minimum 5, larger values provide higher contrast
A=zeros(n);
A_exp=zeros(n);
mod=n; % Select motion mode
for i=1:n
    for j=1:n
        for k=2:mod
            A(i,j)=A(i,j)+V(i,k)*V(j,k)/D(k,k); % ANM three V(i,k)*V(j,k)/D(k,k) summation
            A_exp(i,j)=A_exp(i,j)+(V(i,k)*V(j,k)/D(k,k))*exp(-D(k,k)*tau);
        end
    end
end
T=zeros(n);
for i=1:n
    for j=1:n
        if i~=j
            a=A(j,j)*A(j,j)-A_exp(j,j)*A_exp(j,j);
            b=A(i,i)*A(j,j)*A(j,j);
            c=A(i,j)*A_exp(j,j)*A_exp(i,j);
            d=(A_exp(i,j)*A_exp(i,j)+A(i,j)*A(i,j))*A(j,j);
            e=A_exp(j,j)*A_exp(j,j)*A(i,i);
            f=A(j,j);
            g=A(i,i)*A(j,j)-A(i,j)*A(i,j);
            T(i,j)=0.5*log(a)-0.5*log(b+2.0*c-d-e)-0.5*log(f)+0.5*log(g);
      
       end
    end
end
netTE=T-T';
norm_netTE=netTE./max(max(netTE)); % Normalized net transfer entropy matrix for plotting
norm_difference=sum(norm_netTE,2); % Normalized net transfer entropy for each residue, reflecting emission/reception status
norm_ave = norm_difference./i;
figure(6);
plot(norm_ave);
xlabel('Residue number');
ylabel('net Transfer Entropy');
title('Transfer Entropy');

end