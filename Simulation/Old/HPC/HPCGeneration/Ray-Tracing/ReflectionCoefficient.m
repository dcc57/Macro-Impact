function R = ReflectionCoefficient(LayerId,TempAmp,sinTheta,AlphaRay,v1,v2,rho1,rho2)
%The fluid in which the wave arrives is labeled 1
%   Detailed explanation goes here
    R = 1 - TransmissionCoefficient(LayerId,TempAmp,sinTheta,AlphaRay,v1,v2,rho1,rho2);
    
end

