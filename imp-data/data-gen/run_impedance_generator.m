function run_impedance_generator(test_id)
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

geoinfo = [];
kinfo = [];
impedance_type = 'constkappa';
path_to_ios2d = '../../inverse-obstacle-scattering2d/';
path_to_output_folder = '../data-out/';
verbose = true;

switch test_id
    case 16
        
        % starfish with transmission defaults

%          NOTE: in the AB notation b1 = delta/omega, 
%                     b2 = 1/(rho_r*c_r*sqrt(1+delta^2/omega^2)), and
%                     b3 = 1/(rho_r*(1+delta^2/omega^2))        
        geoinfo.name = 'starfish';
        impedance_type = 'antbar3';
        geoinfo.narm = 5;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        c1 = 0.5;
        c2 = 1.0;
        rho1 = 1.2;
        rho2 = 0.7;
        rhor = rho1/rho2;
        cr = c1/c2;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        delta = sqrt(3)*khmax;
        lamcfs = cell(3,1);
        lamcfs{1} = @(kh) delta/(kh*c2); % delta/omega
        lamcfs{2} = @(kh) 1/(rhor*cr*sqrt(1+(delta/(kh*c2))^2)); 
        lamcfs{3} = @(kh) 1/(rhor*(1+(delta/(kh*c2))^2));
        
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 17
        
        % less complicated plane
        %
        % smooth plane with transmission defaults

%          NOTE: in the AB notation b1 = delta/omega, 
%                     b2 = 1/(rho_r*c_r*sqrt(1+delta^2/omega^2)), and
%                     b3 = 1/(rho_r*(1+delta^2/omega^2))        
        geoinfo.name = 'smooth_plane';
        impedance_type = 'antbar3';
        geoinfo.nterms = 40;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 59;
        c1 = 0.5;
        c2 = 1.0;
        rho1 = 1.2;
        rho2 = 0.7;
        rhor = rho1/rho2;
        cr = c1/c2;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        delta = sqrt(3)*khmax;
        lamcfs = cell(3,1);
        lamcfs{1} = @(kh) delta/(kh*c2); % delta/omega
        lamcfs{2} = @(kh) 1/(rhor*cr*sqrt(1+(delta/(kh*c2))^2)); 
        lamcfs{3} = @(kh) 1/(rhor*(1+(delta/(kh*c2))^2));
        
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 18
        
        % complicated plane
        %
        % smooth plane with transmission defaults

%          NOTE: in the AB notation b1 = delta/omega, 
%                     b2 = 1/(rho_r*c_r*sqrt(1+delta^2/omega^2)), and
%                     b3 = 1/(rho_r*(1+delta^2/omega^2))        
        geoinfo.name = 'smooth_plane';
        impedance_type = 'antbar3';
        geoinfo.nterms = 70;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 79;
        c1 = 0.5;
        c2 = 1.0;
        rho1 = 1.2;
        rho2 = 0.7;
        rhor = rho1/rho2;
        cr = c1/c2;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        delta = sqrt(3)*khmax;
        lamcfs = cell(3,1);
        lamcfs{1} = @(kh) delta/(kh*c2); % delta/omega
        lamcfs{2} = @(kh) 1/(rhor*cr*sqrt(1+(delta/(kh*c2))^2)); 
        lamcfs{3} = @(kh) 1/(rhor*(1+(delta/(kh*c2))^2));
        
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        
    case 19
        
        % random domain
        %
        % with transmission defaults

%          NOTE: in the AB notation b1 = delta/omega, 
%                     b2 = 1/(rho_r*c_r*sqrt(1+delta^2/omega^2)), and
%                     b3 = 1/(rho_r*(1+delta^2/omega^2))        

        rng(1234);
        geoinfo.name = 'random';
        impedance_type = 'antbar3';
        geoinfo.nmode = 7;
        kinfo.k1 = 1;
        kinfo.dk = 0.5;
        kinfo.nk = 39;
        c1 = 0.5;
        c2 = 1.0;
        rho1 = 1.2;
        rho2 = 0.7;
        rhor = rho1/rho2;
        cr = c1/c2;
        khmax = kinfo.k1+kinfo.dk*(kinfo.nk-1);
        delta = sqrt(3)*khmax;
        lamcfs = cell(3,1);
        lamcfs{1} = @(kh) delta/(kh*c2); % delta/omega
        lamcfs{2} = @(kh) 1/(rhor*cr*sqrt(1+(delta/(kh*c2))^2)); 
        lamcfs{3} = @(kh) 1/(rhor*(1+(delta/(kh*c2))^2));
        
        % use the default delta which is just in the 
        % suggested asymptotic regime in Antoine Barucq paper
        

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
