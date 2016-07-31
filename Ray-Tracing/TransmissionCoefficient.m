function T = TransmissionCoefficient(LayerId,sinTheta,AlphaRay,v1,v2,rho1,rho2)
%The fluid in which the wave arrives is labeled 1
%   Detailed explanation goes here
    
    Theta1 = AlphaRay;
    Theta2 = asin(sinTheta);
    Z1 = rho1.* v1;
    Z2 = rho2.* v2;
    Tpre = (2.* Z2.*sec(Theta2))./(Z1.*sec(Theta1)+Z2.*sec(Theta2));
    Thigh = Tpre > 1;
    Tlow = Tpre < 0;
    Tmid = (true(size(Tpre)) & (~ Thigh)) & (~ Tlow);
    T(Thigh) = 1;
    T(Tlow) = 0;
    T(Tmid) = Tpre(Tmid);
    
end

