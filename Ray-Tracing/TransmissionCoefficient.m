function T = TransmissionCoefficient(LayerId,sinTheta,v1,v2,rho1,rho2)
%The fluid in which the wave arrives is labeled 1
%   Detailed explanation goes here
    m = rho2./rho1;
    n = v1./v2;
    Theta = asin(sinTheta);
    T = abs((2.*m.*cos(Theta))./(m.*cos(Theta)+n.*(1-(n.^(-2)).*(sin(Theta)).^(2)).^(1/2)));

end

