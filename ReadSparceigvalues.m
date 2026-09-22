function ReadSparceigvalues(element,tnkpt,nband,fpath,fname)
clc

kpt = ones(tnkpt,3);
eign = zeros(tnkpt,nband);
occ = zeros(tnkpt,nband);

filename = fullfile(fpath,fname);
%filename = strcat(element,'.eigen');
[fid,~] = fopen(filename,'r');
assert(fid~=-1,'Error: Check the element name and the corresponding .eigen file');
textscan(fid,'%s',1,'delimiter','\n');

for i = 1:tnkpt
	textscan(fid,'%s',1,'delimiter','\n');
	textscan(fid,'%s',1,'delimiter','(');
	kpt(i,1) = fscanf(fid,'%f',1);
    textscan(fid,'%s',1,'delimiter',',');
	kpt(i,2) = fscanf(fid,'%f',1);
    textscan(fid,'%s',1,'delimiter',',');
	kpt(i,3) = fscanf(fid,'%f',1);
	textscan(fid,'%s',2,'delimiter','\n');
	for j = 1:nband
		band = fscanf(fid,'%d',1);
		eign(i,band) = fscanf(fid,'%f',1);
        occ(i,band) = fscanf(fid,'%f',1);
		textscan(fid,'%s',1,'delimiter','\n');
	end
end

nu = unique(kpt(:,2));
snu = size(nu,1); % remove factor of 2
kptsplit_ind = ones(snu+1,1);

count = 2;
for i = 2:tnkpt
    if (kpt(i,2) ~= kpt(i-1,2))
        kptsplit_ind(count) = i;
        count = count + 1;
    end
end
kptsplit_ind(end) = tnkpt + 1;

file = strcat('SPARCEigenvalues_',element,'.mat');
file_loc = fullfile(fpath,file);
save(file_loc,'eign','kpt','occ','kptsplit_ind');
    
    
