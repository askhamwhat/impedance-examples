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
path_to_output_folder = '../data-out/';
verbose = true;

switch test_id
    case 1
        % go full defaults...

    case 2
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
    case 3
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 4
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper

    case 5
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2*khmax;
        % half of the default delta so it is just under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 6
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4*khmax;
        % 1/4 of the default delta so it is just under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 7
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8*khmax;
        % 1/8 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 8
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        % 1/16 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 9
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/32*khmax;
        % 1/32 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 10
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/64*khmax;
        % 1/64 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 11
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/128*khmax;
        % 1/128 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 12
        %
        % complicated plane, trying more frequencies
        %
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 13
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2*khmax;
        % half of the default delta so it is just under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 14
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4*khmax;
        % 1/4 of the default delta so it is just under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 15
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8*khmax;
        % 1/8 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 16
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        % 1/16 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 17
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/32*khmax;
        % 1/32 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 18
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/64*khmax;
        % 1/64 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 19
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/128*khmax;
        % 1/128 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 20
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 21
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2*khmax;
        % half of the default delta so it is just under the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 22
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4*khmax;

    case 23
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8*khmax;

    case 24
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        
    case 25
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/32*khmax;
        
    case 26
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/64*khmax;
        
    case 27
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/128*khmax;
        
    case 28
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/256*khmax;
        
    case 29
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/512*khmax;
        
    case 30
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/1024*khmax;
        
    case 31
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2048*khmax;
        
    case 32
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4096*khmax;
        
    case 33
        % smooth plane with transmission defaults
        geoinfo.name = 'starfish';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 19;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8192*khmax;
        
    case 34
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/256*khmax;
        % 1/128 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        

    case 35
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/512*khmax;
        % 1/128 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper
        

    case 36
        % smooth plane with transmission defaults
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 20;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/1024*khmax;
        % 1/128 of the default delta so it is well under the 
        % suggested asymptotic regime in Antoine Barucq paper

    case 37
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at one reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = 0;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        %
        
    case 38
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        %
        
    case 39
        % smooth plane with transmission defaults
        % less dissipation
        % only collect data at one reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = 0;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        %
        
    case 40
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        %
        

    case 41
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at one transmited receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'transmit';
        geoinfo.angle = 0;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        %
        
    case 42
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of transmited receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'transmit';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        %
        
    case 43
        % smooth plane with transmission defaults
        % less dissipation
        % only collect data at one transmited receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'transmit';
        geoinfo.angle = 0;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        %
        
    case 44
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of transmited receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'transmit';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        %
        
    case 45
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2*khmax;
        %
        
    case 46
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4*khmax;
        %
        
    case 47
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8*khmax;
        %
        
    case 48
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/32*khmax;
        %
        
    case 49
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/64*khmax;
        %
        
    case 50
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        transparams.c1 = 0.5;
        transparams.c2 = 1.0;
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/128*khmax;
        %
        
    case 51
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        %
        
    case 52
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/2*khmax;
        %
        
    case 53
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/4*khmax;
        %
        
    case 54
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/8*khmax;
        %
        
    case 55
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/16*khmax;
        %
        
    case 56
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/32*khmax;
        %
        
    case 57
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/64*khmax;
        %
        
    case 58
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
        transparams.rho1 = 0.7;
        transparams.rho2 = 1.2;
        transparams.delta = sqrt(3)/128*khmax;
        %
        
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
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
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

    case 72
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
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

    case 73
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
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

    case 74
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
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

    case 75
        % smooth plane with transmission defaults
        % default dissipaiton
        % only collect data at a section of reflected receptor location per incident
        % wave
        geoinfo.name = 'smooth_plane';
        geoinfo.nterms = 70;
        geoinfo.arrangement = 'reflect';
        geoinfo.angle = pi/4;
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
        
    otherwise
        warning('unknown test, doing nothing');
        return

end

fprintf('running test %d generator ...\n',test_id);
fname = generate_transmission_tensor_data(test_id,transparams,...
    geoinfo,kinfo,path_to_ios2d,path_to_output_folder,verbose);
fprintf('done. output written to: %s\n',fname);
