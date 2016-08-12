function output = QualityFactor(LayerId,DeltaTime,AverageV,Q,Suppression)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
output = Suppression + ((AverageV(LayerId).*DeltaTime)./Q(LayerId));
end