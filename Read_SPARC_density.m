clc; clear; close all;
dens=reader('DPN_00.dens');
dens_avg_z=squeeze(mean(dens,[1 2]));
dens_avg_x = squeeze(mean(dens,[2 3]));
dens_avg_y = squeeze(mean(dens,[1 3]));

Lx = 25.3201285757752;
Ly = 18.1202725592092;
Lz = 28.3720360473642;
nx = 196; ny = 141; nz = 220;
X = linspace(0, Lx, nx);
Y = linspace(0, Ly, ny);
Z = linspace(0, Lz, nz);
tol = 1e-6;

idx_dens = find(abs(dens_avg_z) < tol);
[left_dens, right_dens] = index_finder(idx_dens);

fprintf('density index %d and %d\n',left_dens, right_dens);
fprintf('Z value %.6f and %.6f\n',Z(left_dens), Z(right_dens));

idx_dens = find(abs(dens_avg_x) < tol);
[left_dens, right_dens] = index_finder(idx_dens);

fprintf('density index %d and %d\n',left_dens, right_dens);
fprintf('X value %.6f and %.6f\n',X(left_dens), X(right_dens));

idx_dens = find(abs(dens_avg_y) < tol);
[left_dens, right_dens] = index_finder(idx_dens);

fprintf('density index %d and %d\n',left_dens, right_dens);
fprintf('Y value %.6f and %.6f\n',Y(left_dens), Y(right_dens));

figure;
plot(Z,dens_avg_z,'LineWidth',2);
hold on;
plot(X,dens_avg_x,'LineWidth',2);
plot(Y,dens_avg_y,'LineWidth',2);
xlabel('R'); ylabel('Density');
title('Average rho vs r for DPN');
legend('Averge along Z','Averge along X','Averge along Y');
grid on


function phi = reader(filename)
    % Open file
    fileID = fopen(filename, 'r');

    % Read the second line for cell lengths
    fgetl(fileID);  % Skip first header line
    line2 = fgetl(fileID);  % Second line with cell length info
    tokens = regexp(line2, 'Cell length:\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)', 'tokens');
    
    if isempty(tokens)
        error('Cell length line not in expected format.');
    end
    lengths = str2double(tokens{1});
    Lx = lengths(1);

    % Skip next line, then read nx, ny, nz
    fgetl(fileID);  % Line 3
    nx = sscanf(fgetl(fileID), '%d', 1);  % Line 4
    ny = sscanf(fgetl(fileID), '%d', 1);  % Line 5
    nz = sscanf(fgetl(fileID), '%d', 1);  % Line 6

    % Skip the next few metadata lines (assumed until line 10)
    while true
        line = fgetl(fileID);
        if isnumeric(line) || isempty(line)
            continue;
        end
        vals = sscanf(line, '%f');
        if length(vals) ~= 5
            break;
        end
    end

    % Go back to start of potential data line
    fseek(fileID, -numel(line) - 1, 'cof');

    % Read potential data
    N = nx * ny * nz;
    data = fscanf(fileID, '%f');
    data = data(1:N);
    % Reshape to 3D grid
    phi = reshape(data, [nz, ny, nx]);
    phi = permute(phi, [3, 2, 1]);  % Convert to (nx, ny, nz)
end

function [left_index, right_index] = index_finder(idx_dens)
   idx = idx_dens(:);
    breaks = find(diff(idx) > 1);
    if isempty(breaks)
        left_index  = idx(end);
        right_index = NaN;
    else
        left_index  = idx(breaks(1));
        right_index = idx(breaks(1) + 1);
    end
end
