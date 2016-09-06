function output = T(x,y,z,theta)
%Rotation about the vector (ux,uy,uz) by the angle theta.
%   Detailed explanation goes here
    mag = sqrt(x^2+y^2+z^2);
    ux = x/mag;
    uy = y/mag;
    uz = z/mag;
    t = 1 - cos(theta);
    C = cos(theta);
    S = sin(theta);
    output = [t*(ux^2)+C,t*ux*uy-S*uz,t*ux*uz+S*uy;t*ux*uy+S*uz,t*(uy^2)+C,t*uy*uz-S*ux;t*ux*uz-S*uy,t*uy*uz+S*ux,t*(uz^2)+C];
end

