clear

format long

R = 1.737 * 10^6; % m
vsurf = 1000; %m/s
Res = 1000; %Number of subdivisions per dT

load /Users/David/Macro-Impact/Simulation/Current/Cart2Hp/hppoints1.mat hppoints
load /Users/David/Macro-Impact/Simulation/Current/DataGeneration/Mat1.mat outputmat
load /Users/David/Macro-Impact/Simulation/Current/DataGeneration/Num1.mat outputnum

hppoints = double(hppoints);
outputmat = double(outputmat);
outputnum = double(outputnum);

Data = cat(2,transpose(hppoints),outputmat(:,4),outputmat(:,5),outputmat(:,6),outputmat(:,8),outputmat(:,9)); % t,rho,REF,kappa,vp

m = outputnum(1);
N = outputnum(2);
L = outputnum(3);

clear('hppoints','outputmat','outputnum')
A = 4 * pi * R^2 / N;
dT = sqrt(A / pi) / vsurf;

Data = sortrows(Data,1);
C = 1;
Disp = zeros(N,1);

%NEW DATA ANALYSIS (MUCH FASTER)

%%{
Tmax = max(Data(:,2));
Tmin = min(Data(:,2));

Bins = int64(Res * ceil((Tmax-Tmin)/dT));

T = zeros(length(Data(:,2)),1);

for i = 1 : Bins
    Id = (Data(:,2) > Tmin + (i - 1) * dT / Res) & (Data(:,2) <= Tmin + i * dT / Res);
    T(Id,1) = i;
end

Data(:,2) = T;
Data = sortrows(Data,2);

Points = Data(:,1);
T = Data(:,2);

TempDisp = zeros(N,Bins);
TempDispMax = zeros(N,1);

AllTempDisp = (Data(:,3) .* Data(:,4) .* Data(:,5)) ./ (Data(:,6).^3);
SumTempDisp = zeros(N,Bins);
parfor j = 0 : N - 1
    for i = 1 : Bins
        Id = (T == i) & (Points == j);
        SumTempDisp(j + 1,i) = sum(AllTempDisp(Id));
    end
end

Numbers = transpose(1 : Bins);

for j = 0 : N - 1;
    for i = 1 : Bins
        Id = Numbers >= i - 1 & Numbers < Res + i;
        TempDisp(j + 1,i) = sum(SumTempDisp(j + 1,Id));
    end
    TempDispMax(j + 1) = max(TempDisp(j + 1,:));
end

PreDisp = ((TempDispMax * 9 * L) / (2 * pi * N * m * A)).^(1/2);

%}

%OLD DATA ANALYSIS (FOR COMPARISON)

%{
tic

for i = 0 : N - 1
    I = nnz(Data(:,1) == i);
    TempData = Data(C : C + I - 1,2 : 6);
    C = C + nnz(Data(:,1) == i);
    for j = 1 : I
        LID = abs(TempData(j,1) - TempData(:,1)) < dT;
        Disp(i+1,1) = max(Disp(i+1,1),sum((TempData(LID,2) .* TempData(LID,3) .* TempData(LID,4)) ./ (TempData(LID,5).^3)));
    end
end

toc

PreDisp = ((Disp * 9 * L) / (2 * pi * N * m * A)).^(1/2);
%}

save PreDisp1.mat PreDisp;