% --------------------------------------------------------
% plot_tdkp_gain
%
% plot gain data produced by tdkp / clc
% --------------------------------------------------------
function [omega,Y] = plot_tdkp_gain(clc_mode, densities_file, source_dir, column, linestyle, minus)


    if nargin == 0
        disp('usage: plot_tdkp_gain(clc_mode, densities_file, source_dir, column, linestyle, minus)');
        return
    end


    if nargin <= 4
        linestyle = 'r-';
    end
    if nargin <= 5
        minus = 1.0
    end

    [omega, Y] = load_tdkp_gain(clc_mode, column, densities_file, source_dir);

    Y = Y * minus;

    plot(omega,Y, linestyle);







% ------------------------------------
% function load_tdkp_gain(densities_file,sourcedir)
%
% loads gain data from multiple files
% ------------------------------------
function [w, Y] = load_tdkp_gain(clc_mode, column, densities_file,sourcedir)

    dens   = load_tdkp_data(densities_file);
    Y = [];
    for ii = 1:size(dens,1)
        if clc_mode == 1
            dt = load_tdkp_data(sprintf('%s/clc_optics_n%g_p%g.dat',sourcedir,dens(ii,1), dens(ii,2)));
        else
            dt = load_tdkp_data(sprintf('%s/optics_n%g_p%g.dat',sourcedir,dens(ii,1), dens(ii,2)));
        end
        Y(:,ii) = dt(:,column);
    end
    w = dt(:,1);


% ------------------------------------
% function load_tdkp_data(filename)
%
% similar to load -ascii except that
% comments (lines starting with #)
% are ignored
% ------------------------------------
function X = load_tdkp_data(filename)

    fid = fopen(filename);
    X   = [];
    ii  = 1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if tline(1) ~= '#'
            t = sscanf(tline, '%g');
            if ii > 1 && length(t) ~= size(X,2)
                error('unequallength',sprintf('number of values in line %d is %d while i expected %d', ii, length(t), size(X,2)));
            end
            X(ii,:) = t(:);
            ii = ii + 1;
        end
    end
    fclose(fid);
    return
