function T = TransmissionCoefficient(LayerId,TempAmp,sinTheta,AlphaRay,v1,v2,rho1,rho2)
%The fluid in which the wave arrives is labeled 1
%   Detailed explanation goes here
    Theta1 = zeros(size(AlphaRay));
    Id = AlphaRay < pi / 2;
    Idp = AlphaRay >= pi / 2;
    Theta1(Id) = AlphaRay(Id);
    Theta1(Idp) = pi - AlphaRay(Idp);
    Theta2 = asin(sinTheta);
    Z1 = rho1.* v1;
    Z2 = rho2.* v2;
    Tpre = (2.* Z2.*sec(Theta2))./(Z1.*sec(Theta1)+Z2.*sec(Theta2));
    Thigh = Tpre > 1;
    Tlow = Tpre < 0;
    Tmid = (true(size(Tpre)) & (~ Thigh)) & (~ Tlow);
    T(Thigh) = TempAmp(Thigh);
    T(Tlow) = 0;
    T(Tmid) = Tpre(Tmid).*TempAmp(Tmid);
    
end

