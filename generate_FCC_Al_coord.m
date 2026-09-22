clc; clear;
function frac_coords = generate_fcc_frac(nx, ny, nz)
    % FCC conventional fractional basis
    base = [
        0.0 0.0 0.0;
        0.5 0.5 0.0;
        0.0 0.5 0.5;
        0.5 0.0 0.5
    ];
    frac_coords = [];
    % loop over supercell replicas
    for i = 0:nx-1
        for j = 0:ny-1
            for k = 0:nz-1
                shift = [i j k];
                new_atoms = (base + shift) ./ [nx ny nz];
                frac_coords = [frac_coords; new_atoms];
            end
        end
    end
    % wrap into [0,1)
    frac_coords = mod(frac_coords,1);
end

generate_fcc_frac(8,8,8)
