clear

format long

R = 1.737 * 10^6; % m
vsurf = 1000; %m/s

load hppoints10.mat hppoints
load PointsTestMat10.mat outputmat
load PointsTestNum10.mat outputnum

hppoints = double(hppoints);
outputmat = double(outputmat);
outputnum = double(outputnum);

Data = cat(2,transpose(hppoints),outputmat(:,4),outputmat(:,5),outputmat(:,6),outputmat(:,8),outputmat(:,9)); % t,rho,REF,kappa,vp

m = outputnum(1);
n = outputnum(2);
L = outputnum(3);

clear('hppoints','outputmat','outputnum')
N = n;
A = 4 * pi * R^2 / N;
dT = sqrt(A / pi) / vsurf;

Data = sortrows(Data,1);
C = 1;
Disp = zeros(N,1);

for i = 0 : N - 1
    I = nnz(Data(:,1) == i);
    TempData = Data(C : C + I - 1,2 : 6);
    C = C + nnz(Data(:,1) == i);
    for j = 1 : I
        LID = abs(TempData(j,1) - TempData(:,1)) < dT & abs(TempData(j,1) - TempData(:,1)) > 0;
        Disp(i+1,1) = max(Disp(i+1,1),sum((TempData(LID,2) .* TempData(LID,3) .* TempData(LID,4)) ./ (TempData(LID,5).^3)));
    end
end

PreDisp = ((Disp * 9 * L) / (2 * pi * n * m * A)).^(1/2);

save PreDisp10.mat PreDisp;