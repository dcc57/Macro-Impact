clear

format long
SampleSize = 17;
Threshold = 3E-10;

fmax = 20;
vX = 2.7E5;
sX = 1E-11;
p0 = 3.0E8;

C = zeros(SampleSize,1);
D = zeros(SampleSize,1);

for i = 1 : SampleSize
    A = strcat('PreDisp',int2str(i));
    B = strcat('PointsTestNum',int2str(i));
    load(A,'PreDisp');
    load(B,'outputnum');
    m = outputnum(1);
    n = outputnum(2);
    L = outputnum(3);
    D(i) = outputnum(4);
    Disp = PreDisp .* ((p0 ^ (-1)) * (fmax ^ (1/2)) * sX * (vX ^ 2));
    C(i) = nnz(Disp > Threshold);
end
F = (sum(C.*D)/sum(D))/n