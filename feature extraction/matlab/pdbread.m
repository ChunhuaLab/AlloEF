function [posall,resall,bf]=pdbread(filename)
fclose('all');
max=35000;
posall=[];
resall=[];
bf=[];
ft = fopen(filename, 'r');
if ft == -1
    error('Cannot open file: %s', filename);
end
for i=1:max
    pl=fgetl(ft);
    if pl(1)=='E'& pl(2)=='N'& pl(3)=='D'
        break;
    elseif ((pl(1:4)=='ATOM')&((pl(14)=='C')&(pl(15)=='A')))|((pl(1:6)=='HETATM')&((pl(14)=='C')&(pl(15)=='A')))
        rn=pl(18:20);
        if rn=='GLY'
            resname='G';
        elseif rn=='ALA'
            resname='A';
        elseif rn=='VAL'
            resname='V';
        elseif rn=='LEU'
            resname='L';
        elseif rn=='ILE'
            resname='I';
        elseif rn=='PRO'
            resname='P';
        elseif rn=='PHE'
            resname='F';
        elseif rn=='TRP'
            resname='W';
        elseif rn=='TYR'
            resname='Y';
        elseif rn=='SER'
            resname='S';
        elseif rn=='THR'
            resname='T';
        elseif rn=='CYS'
            resname='C';
        elseif rn=='MET'
            resname='M';
        elseif rn=='ASN'
            resname='N';
        elseif rn=='GLN'
            resname='Q';
        elseif rn=='ASP'
            resname='D';
        elseif rn=='GLU'
            resname='E';
        elseif rn=='HIS'
            resname='H';
        elseif rn=='LYS'
            resname='K';
        elseif rn=='ARG'
            resname='R';
        else
            resname='U';
        end
        resall=[resall resname];
        pos=str2num(pl(30:54));
        posall=[posall;pos];
        bfa=str2num(pl(61:66));
        % bf=[bf;bfa];
        if bfa < 1.00 % B-factor < 1 in homology modeling part has issues, better to use original pdb file to find max PCC parameters, use complete homology model structure for network building
            bf=[bf;bfa*100];
        else
            bf=[bf;bfa];
        end
    end
end
fclose('all');