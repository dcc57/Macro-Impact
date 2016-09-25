clear

format long

R = 1.737 * 10^6; % m
vsurf = 1000; %m/s

load hppoints1Test.mat hppoints
load PointsTestMat1Test.mat outputmat
load PointsTestNum1Test.mat outputnum

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
IterDataLength = 0;
Out = zeros(N,1);

SumDisp = zeros(N,1);

Disp = (Data(:,3) .* Data(:,4) .* Data(:,5)) ./ (Data(:,6).^3);
Time = Data(:,2);

for i = 1 : N
    IterDataLength = max(IterDataLength,length(Data(Data(:,1) == i)));
end
IterDisp = zeros(IterDataLength,N);
IterTime = zeros(IterDataLength,N);
for i = 1 : N
    IterDisp(:,i) = cat(1,Disp(Data(:,1) == i),zeros(IterDataLength-length(Disp(Data(:,1) == i)),1));
    IterTime(:,i) = cat(1,Data((Data(:,1) == i),2),zeros(IterDataLength-length(Disp(Data(:,1) == i)),1));
end
for i = 1 : N
    Comb = combnk(1 : IterDataLength,2);
    LID = abs(IterTime(Comb(:,1)) - IterTime(Comb(:,2))) < dT & abs(IterTime(Comb(:,1)) - IterTime(Comb(:,2))) > 0;
    %CombIterDisp = combnk(IterDisp(:,i),2);
    for ii = 1 : IterDataLength
        if LID(ii) == 1
            IterOut(ii) = IterDisp(comb(ii,1));
        end
        PixOut(i) = max(IterOut)
    end
end

%parfor i = 1 : N
%    sum(combnk(IterData(LID),2),1)
%end
    
%(TempData(LID,2) .* TempData(LID,3) .* TempData(LID,4)) ./ (TempData(LID,5).^3))

PreDisp = ((Disp * 9 * L) / (2 * pi * n * m * A)).^(1/2);

save PreDisp17.mat PreDisp;