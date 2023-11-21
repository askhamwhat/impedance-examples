
function [delta,rhor,cr] = convert_lamabv_beta_to_phys(betas,omega)
%
% sqrt(1-1i*delta/omega)/(rhor*cr*sqrt(1+delta^2/omega^2)) = ...
%              betas(2)*sqrt(1-1i*betas(1))
% (1-1i*delta/omega)/(rhor*(1+delta^2/omega^2)) = betas(3)*(1-1i*betas(1))
%
% delta = omega*betas(1)
% rhor = 1/(betas(3)*(1+delta^2/omega^2))
% cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2)) 

delta = omega*betas(1);
rhor = 1/(betas(3)*(1+betas(1)^2));
cr = 1/(rhor*betas(2)*sqrt(1+delta^2/omega^2));

end