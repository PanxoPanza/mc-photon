function [R, T, r, t] = Fresnel(theta,n1,n2,pol)
% Reflection and transmission based on Fresnel coefficients
% Input:
%   TT: angle of incidece (in degrees)
%   nn1: Refractive index of upper material
%   nn2: Refractive index of lower material
% Output:
%   R: Refleted power
%   T: Transmitted power
%   r: reflection coefficient
%   t: reflection coefficient

if (theta > 90) 
    disp('Error incidence angle must be less than 90 deg');
    return
end

% [Tt,nn1] = meshgrid(theta,n1);
% [~,nn2] = meshgrid(theta,n2);

[nn1, Tt] = meshgrid(n1,theta);
nn2       = meshgrid(n2,theta);

Tt = pi*Tt/180;
sinTi = sin(Tt);
cosTi = cos(Tt);
sinTt = nn1.*sinTi./nn2;
cosTt = sqrt(1 - sinTt.^2);
if (~strcmp(pol,'p') && ~strcmp(pol,'s'))
    disp('Polarization not defined. Use p or s');
    return
else 
    if strcmp(pol,'p')
        r = (nn1.*cosTt - nn2.*cosTi)./ (nn1.*cosTt + nn2.*cosTi);
        t = (2*nn1.*cosTi)           ./ (nn1.*cosTt + nn2.*cosTi);
        R = r.*conj(r);
        T = real(conj(nn2).*cosTt)./real(conj(nn1).*cosTi).*t.*conj(t);
        
    elseif strcmp(pol,'s')
        r = (nn1.*cosTi - nn2.*cosTt)./ (nn1.*cosTi + nn2.*cosTt);
        t = (2*nn1.*cosTi)           ./ (nn1.*cosTi + nn2.*cosTt);
        R = r.*conj(r);
        T = real((nn2).*cosTt)./real((nn1).*cosTi).*t.*conj(t);
    end
end
R = abs(R); T = abs(T);
end