function output = QualityFactor(LayerId,DeltaTime,vBoundaries,Q,Suppression)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
output = Suppression + ((vBoundaries(LayerId).*DeltaTime)./Q(LayerId));
end