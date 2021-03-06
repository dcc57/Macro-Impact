clear

format long

Save = 1; %1 = save displacement arrays, 0 = do not save
SampleSize = 1;
Threshold = 1* 3E-10;
p = log(sqrt(2)); %Detection Likelihood Scaling

fmax = 20;
vX = 2.5E5;
sX = 1E-11;
p0 = 1.0E8;

C = zeros(SampleSize,1);
Cweight = zeros(SampleSize,1);
D = zeros(SampleSize,1);

for i = 1 : SampleSize
    A = strcat('/Users/David/Macro-Impact/Simulation/Current/DataAnalysis/PreDisp',int2str(i));
    B = strcat('/Users/David/Macro-Impact/Simulation/Current/DataGeneration/Num',int2str(i));
    load(A,'PreDisp');
    load(B,'outputnum');
    m = outputnum(1);
    n = outputnum(2);
    L = outputnum(3);
    D(i) = outputnum(4);
    Disp = PreDisp .* ((p0 ^ (-1)) * (fmax ^ (1/2)) * sX * (vX ^ 2));
    WeightedDisp = max(1 - exp(p.*(-(Disp-Threshold)./Threshold)),0);
    C(i) = nnz(Disp> Threshold);
    Cweight(i) = sum(WeightedDisp);
    
    %SAVE DISPLACEMENTS TO MAKE PICTURES
    
    if Save == 1;
        S = strcat('Disp',int2str(i),'.mat');
        save(S,'Disp');
    end
end
F = (sum(C.*D)/sum(D))
Fweight = (sum(Cweight.*D)/sum(D))