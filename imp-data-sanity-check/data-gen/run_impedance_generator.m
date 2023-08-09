%RUN_IMPEDANCE_GENERATOR
%
% this file has the parameters for various test runs 
% 
% Note: there isn't much relation between test number and test difficulty
% because of debugging
%
% these are for constant + kappa type models. see the function 
% /src/matlab/constkappa_models_convert for details 
%
% constkappa model has two complex parameters 
% antbar2 has two real parameters (assumes rho_r = 1)
% antbar3 has three real parameters
clearvars; 

test_id = 13;

geoinfo = [];
kinfo = [];
impedance_type = 'constkappa';
path_to_ios2d = '../../inverse-obstacle-scattering2d/';
path_to_output_folder = '../data-out/';
verbose = true;

switch test_id
    case 1
        lamcfs = [0.5-0.3*1i; 0.2 + 0.3*1i];
    case 2
        geoinfo.name = 'starfish';
        geoinfo.narm = 0;
        geoinfo.rad = 1.1;
        lamcfs = [1.1+0.9*1i,0.3+0.2*1i];
    case 3
        geoinfo.name = 'starfish';
        geoinfo.narm = 0;
        geoinfo.rad = 1.5;
        lamcfs = [1.1-0.9*1i,0];
    case 4
        geoinfo.name = 'starfish';
        geoinfo.narm = 2;
        geoinfo.rad = 1.0;
        lamcfs = [1.1-0.9*1i,0.1+0.2*1i];
    case 5
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 10;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5-0.3*1i; 0.2 + 0.3*1i];
    case 6
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5-0.3*1i; 0.2 + 0.3*1i];
    case 7
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 30;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5-0.3*1i; 0.2 + 0.3*1i];
    case 8
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 50;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5-0.3*1i; 0.2 + 0.3*1i];
    case 9
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5; 0.2];
    case 10
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 30;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5; 0.2];
    case 11
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 50;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [0.5; 0.2];
    case 12
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        impedance_type = 'antbar2';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        lamcfs = [3.0; 1.0];
    case 13
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        impedance_type = 'antbar2';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        delta = 20; % while this seems large it gets attenuated by omega
        c2 = 2.1;
        c1 = 1.2;
        cr = c1/c2;
        % specify omega-dependency by kh = k2 = omega/c2 dependency
        lamcfs = cell(2,1);
        lamcfs{1} = @(kh) delta/(kh*c2); % delta/omega
        lamcfs{2} = @(kh) 1/cr; 
        
    otherwise
        warning('unknown test, doing nothing');
        lamcfs = [];

end

if ~isempty(lamcfs)
    fprintf('running test %d ...\n',test_id);
    fname = generate_impedance_tensor_data(test_id,lamcfs,...
        impedance_type,geoinfo,kinfo,path_to_ios2d, ...
        path_to_output_folder,verbose);
    fprintf('done. output written to: %s\n',fname);
end
