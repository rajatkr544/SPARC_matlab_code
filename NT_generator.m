function [xcyc,xcyc_red,xcart] = NT_generator(a,del,n,m,twist_ext,ncell_c,ncell_z,vacc,bond_typ,class,isPlot)
%%

%==========================================================================
%                         Purpose
%==========================================================================
% Generate atom positions for achiral and chiral nanotubes of different classes:
% (1) in helical coord. for the unit cell with cyclic and helical symmetry 
% (2) in cartesian coord. for the unit cell with only periodic symmetry
% along axis
%%

%==========================================================================
%                       Details
%==========================================================================
% (1) Nanotubes are generally identified by a pair of indices (n,m)
% (2) Achiral nanotubes are of two types, namely, a) Armchair b) Zig-zag
% (3) Armchair nanotubes are represented by (n,n) (hexagonal) and by (0,m)
% (rectangular)
% (4) Zig-zag nanptubes are represented by (n,0) notation
% (5) 4-atom cell typ for hexagonal class is available only for achiral tubes
% (6) Twist (rad/length) can be given for both achiral and chiral tubes
% (7) Screw translation is represented by {psi|tau}
%%

%==========================================================================
%                      Input to the code
%==========================================================================
% (1a) a       :       X-X bond length
% (1b) del     :       out of plane projection of the atom
% (2a) (n,m)   :       indices representing the type of nanotube
% (2b)twist_ext:       gives the applied twist per unit length
% (3) ncell_c  :       number of cells in the cyclic direction (only for
% Cyclix code)
% (4) ncell_z  :       number of cells in the axial direction (only for
% Cyclix code)
% (5) vacc     :       vacuum needed in the radial direction
% (6) bond_typ :       type of bond
%                      1 = linear, 2 = arc
% (7) class    :       gives the class of the nanotube (X2,XY,XY2)
% (8) isPlot   :       plot the nanotube
%%

%==========================================================================
%                      Output from the code
%==========================================================================
% xcyc         :   helical coord. of atoms in ncell_c x ncell_z cells
% xcyc_red     :   reduced coordinates corresponding to xcyc coordinates
% xcart        :   cartesian coord. of atoms in the 1D periodic cell
%%
clc
format long

% Unit conversion
bohr2ang = 0.529177210903;
ang2bohr = 1/bohr2ang;

%==========================================================================
%                        Nanotube characterization
%==========================================================================
A = sqrt(3)*a;
a1 = [sqrt(3)/2 1/2]*A; a2 = [sqrt(3)/2 -1/2]*A;
d = gcd(n,m);
dR = gcd(2*n+m, 2*m+n);
Ch = A * sqrt(n^2 + m^2 + n*m); % Circumference
t1 = (2*m+n)/dR; t2 = -(2*n+m)/dR; 

% (1) Diameter of the tube
if (d == 1 && bond_typ == 1)
    error('Bond_typ = 1 not allowed for systems with no cyclic symmetry');
elseif (d < 10 && bond_typ == 1)
    fprintf('Warning: Bond_typ = 1 may cause a lot of deviation from actual shape\n');
end
if bond_typ == 1
    dia = 2*(Ch/2/d)/sin(pi/d);
else
    dia = Ch/pi;
end

fprintf('Radius of the tube is %.15f Bohr\n',dia/2);
fprintf('Angle subtended by the unit cell %.15f\n',2*pi/d);

% (2) Height of the periodic cell
T = sqrt(3) * Ch /dR;
fprintf('Height of the periodic cell (from helical estimate) is %.15f Bohr\n',T);

% (3) Chiral angle
Ch_theta = acos((2*n+m)/(2*sqrt(n^2+m^2+n*m)));
fprintf('Chiral angle is %.15f radians(%.15f degree)\n',Ch_theta,Ch_theta*180/pi);

% (4) Symmetry operations for 2 atom cells (screw translation)
N_units = 2*(n^2+m^2+n*m)/dR
tau = d*T/N_units;

% Find q under conditions: 
%   a) pm-qn = d, p,q are integers
%   b) -mM_1/N < q < m - mM_1/N

flag = 0; 
if m == 0
    M = 1;
    flag = 1;
else
    M_1 = (2*n+m)*d/m/dR;
    lowerb = floor(-m*M_1/N_units + 1); % lower bound
    upperb = ceil(m - m*M_1/N_units - 1); % upper bound
    for q = lowerb:upperb
        M = M_1 + N_units*q/m;
        p = (d + q*n)/m;
        if (round(M) == M && round(p) == p)
            flag = 1;
            break;
        end
    end
end
M
if flag == 1
    psi = 2*pi*M/N_units;
else
    error('Failed to get M\n');
end
fprintf('Screw translation corresponding to 2 atom cell is {%f|%f}\n',psi,tau);

% (5) Twist (per unit length) to be given to the code for running 2-atom cell simulation
twist_int = psi/tau;
fprintf('Twist of the tube corresponding to 2-atom unit cell %f rad/bohr\n',twist_int);
twist = twist_int+twist_ext;
fprintf('Total twist is %.15f rad/bohr\n',twist);
%%

%==========================================================================
%               Atom positions in rectangular flat sheet
%==========================================================================
cell = [Ch/d tau];
if strcmp(class, 'X2') || strcmp(class,'XY')
	[pos_Ch, pos_T] = findAtom(n,m,a1,a2,t1,t2,Ch,T,cell,N_units);
    x_uc = [0 0 0;pos_Ch/cell(1) pos_T/cell(2) del];
    natom_u = 2;
elseif strcmp(class,'XY2')
    [pos_Ch, pos_T] = findAtom(n,m,a1,a2,t1,t2,Ch,T,cell,N_units);
    x_uc = [0 0 0;pos_Ch/cell(1) pos_T/cell(2) -del;pos_Ch/cell(1) pos_T/cell(2) del];
    natom_u = 3;
end
%%

%==========================================================================
%           Atom positions in cylindrical coordinates (and fractional)
%==========================================================================
rad = dia/2;
if(ncell_c > d)
    error('Number of cells in the cyclic direction (ncell_c) cannot be greater than %d\n',d);
end

xcyc = zeros(natom_u*ncell_c*ncell_z,3);
xcyc(1:natom_u,:) = [rad+x_uc(:,3) x_uc(:,1)*(2*pi/d) x_uc(:,2)*cell(2)]
gen_c = [0 2*pi/d 0];
gen_z = [0 psi tau];
for i = 1:ncell_c-1
    xcyc(natom_u*i+1:natom_u*(i+1),:) = xcyc(1:natom_u,:) + gen_c*i;
end
for i = 1:ncell_z-1
    xcyc(natom_u*ncell_c*i+1:natom_u*ncell_c*(i+1),:) = xcyc(1:natom_u*ncell_c,:) + gen_z*i;
end

% xcyc(:,2) = mod(xcyc(:,2),2*pi);

xmin_at = min(xcyc(:,1));
xmax_at = max(xcyc(:,1));
if(xmin_at - vacc <= 0.1)
    error('vacuum asked too large for the given nanotube!');
end

xcyc_red = xcyc./repmat([xmax_at-xmin_at+2*vacc 2*pi*ncell_c/d cell(2)*ncell_z],size(xcyc,1),1);

fprintf('Fractional coordinates in cylindrical coordinates\n');
disp(xcyc_red);

% Find coordinates corresponding to helical coordinate system
xcyc(:,2) = xcyc(:,2) - twist_int*xcyc(:,3);
xcyc(:,2) = mod(xcyc(:,2),2*pi*ncell_c/d);
xcyc_red(:,2) = xcyc(:,2)/(2*pi*ncell_c/d);

%%

%===============================================================================
%      Atom positions in cartesian coordinates for the 1D periodic cell (only fundamental)        
%===============================================================================
twist = twist_ext+twist_int;
    
% flag = 0;
% for i = 1:10000*d
%     qq = 2*pi*i/d/twist/tau;
%     if abs(round(qq)-qq) < 1e-3
%         ncell_z_cart = round(qq);
%         flag = 1;
%         break;
%     end
% end

ncell_z_cart = N_units/d;

if flag == 1
    fprintf('ncell_z_cart = %d\n',ncell_z_cart);
else
    error('No periodic cell could be found in %f rotation\n',2*pi*10000);
end
xcart = zeros(natom_u*d*ncell_z_cart,3);
for ii = 1:ncell_z_cart
    indx = natom_u*d*(ii-1);
    z0 = (ii-1)*cell(2);
    for i = 1:d
        thetatemp = xcyc(1:natom_u,2) + (i-1)*2*pi/d;
        thetatemp = thetatemp + twist * (z0 + xcyc(1:natom_u,3));
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,1) = xcyc(1:natom_u,1).*cos(thetatemp(:));%xcyc(1:natom_u,1);
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,2) = xcyc(1:natom_u,1).*sin(thetatemp(:));%thetatemp*5;
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,3) = z0 + xcyc(1:natom_u,3);
    end
end
%%

%========================================================================================
%                        Visualization
%========================================================================================
if isPlot == 1 % tube visualization
    %xcart = xcart*bohr2ang;
    %data1.X = (xcart(:,1))'; data1.Y = (xcart(:,2))'; data1.Z = (xcart(:,3))';
    %mat2pdb(data1)
    %molviewer('mat2PDB.pdb')
    %xcart = xcart*ang2bohr;
    sz = size(xcart,1);
    xcart1 = [xcart;xcart;xcart;xcart;xcart;xcart;xcart;xcart;xcart];
    for kk = 1:8
        xcart1(kk*sz+1:(kk+1)*sz,3) = xcart(:,3) + kk*T;
    end
    xcart = xcart1;


    fid = fopen('Nanotube.xyz', 'w');
    fprintf(fid, '%d\n',size(xcart,1));
    fprintf(fid, 'Nanotube\n');
    for i = 1:2:size(xcart,1)
        fprintf(fid,'C %f %f %f\n',xcart(i,:));
        fprintf(fid,'O %f %f %f\n',xcart(i+1,:));
    end
    fclose(fid);
end    
%%

%=============================================================================================
%                       Input for Cyclix SPARC
%=============================================================================================
fprintf('\n\n********************** SPARC INPUT BEGINS (in helical coordinate system)*************\n\n');
fprintf('CELL: %.15f %.15f %.15f\n',xmax_at-xmin_at+2*vacc,2*pi*ncell_c/d,cell(2)*ncell_z);
fprintf('TWIST_ANGLE: %.15f \n',twist);
fprintf('KPT: %d %d %d\n',1,d,ncell_z_cart);
fprintf('[\b"Warning:Change the kpt in the z-direction appropriately"\n]\b');
fprintf('COORD_FRAC:\n')
disp(xcyc_red);
fprintf('************************** SPARC INPUT ENDS ************************************\n\n');
%%

end




function [pos_Ch, pos_T] = findAtom(n,m,a1,a2,t1,t2,Ch,T,cell,N_units)
    flag = 0;
    for p = 0:N_units
        for q = -N_units:N_units
            pos_Ch = dot((p+1/3)*a1 + (q+1/3)*a2,n*a1+m*a2)/Ch;
            pos_T = dot((p+1/3)*a1 + (q+1/3)*a2,t1*a1+t2*a2)/T;
            if pos_Ch + 10*eps >= 0 && pos_Ch < cell(1) && pos_T + 10*eps >= 0 && pos_T < cell(2)
                flag = 1;
                break;
            end
        end
        if(flag == 1)
            break;
        end
    end
    p
    q
    
    if flag == 0
        error('Could not find the atom positions in the 2D sheet\n');
    end
end