function run_transmission_generator(test_id)
%RUN_TRANSMISSION_GENERATOR
%
% this file has the parameters for various test runs 
% 
% Note: there isn't much relation between test number and test difficulty
% because of debugging
%


transparams = [];
geoinfo = [];
kinfo = [];
path_to_ios2d = '../../inverse-obstacle-scattering2d/';
path_to_chunkie = '../../chunkie/';
path_to_output_folder = '../data-out/';
verbose = true;

switch test_id
    case 59
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax;
        
    case 60
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/4;
        
    case 61
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/16;
        
    case 62
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/64;
        
    case 63
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/256;
        
    case 64
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax*4;
        
    case 65
        % smooth plane with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax*16;

    case 66
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/8;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax;
        %

    case 67
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/8;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/4;
        %

    case 68
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/8;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/16;
        %

    case 69
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/8;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/64;
        %

    case 70
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/8;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/256;
        %

    
    case 71
        % L shaped polygon with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'polygon';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax;
        
    case 72
        % L shaped polygon with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'polygon';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/4;
        
    case 73
        % L shaped polygon with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'polygon';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/16;
        
    case 74
        % L shaped polygon with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'polygon';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/64;
        
    case 75
        % L shaped polygon with transmission defaults
        % more physical parameters for sound waves 
        % (speed = bulk modulus/density)
        %

        geoinfo.name = 'polygon';
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 1.2;
        transparams.rho2 = 0.7;
        transparams.delta = sqrt(3)*khmax/256;
        
        
    otherwise
        warning('unknown test, doing nothing');
        return

end

fprintf('running test %d generator ...\n',test_id);
if isfield(geoinfo,'name') && strcmpi(geoinfo.name,'polygon')
    fprintf('this is a chunkie test (polygon)...\n');
    fname = chunkie_corners_generate_transmission_tensor_data(test_id,...
        transparams,geoinfo,kinfo,path_to_chunkie,...
        path_to_output_folder,verbose);
else
    fname = generate_transmission_tensor_data(test_id,transparams,...
        geoinfo,kinfo,path_to_ios2d,path_to_output_folder,verbose);
end
fprintf('done. output written to: %s\n',fname);
