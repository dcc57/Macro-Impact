function output = QualityFactor(LayerId,DeltaTime,Q,Suppression)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
output = Suppression + (DeltaTime./Q(LayerId));
end

