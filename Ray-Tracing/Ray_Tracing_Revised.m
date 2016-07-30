clear

format long

a = [11.2633,10.6835,15.4938].^(1/2); %first entry corresponds to the core
b = [1.55946*10^(-7),2.11853*10^(-7),1.40336*10^(-7)].^(1/2); %first entry corresponds to the core
WeaknessOfNeglectedRays = 10^(-1); %input('What fraction of the initial ray energy must a ray contain to not be thrown out? ');

D = 2000;
Radii = [0,1221,3480,6371];
ACO = 2;

Layers = length(Radii);

vBoundaries = zeros(Layers-2,2); %v(:,1) corresponds to the inner velocity on the :th layer, and v(:,2) corresponds to the outer velocity on the :th layer
ReflectionCoefficients = zeros(Layers - 1,2); %these correspond to 4 functions of the angle of incidence
TransmissionCoefficients = zeros(Layers - 1,2); %these correspond to 4 functions of the angle of incidence
for i = 1 : Layers - 2
    vBoundaries(i,1) = (a(i)^2)-(b(i)^2)*(Radii(i+1))^2; %inner
    vBoundaries(i,2) = (a(i+1)^2)-(b(i+1)^2)*(Radii(i+1))^2; %outer
end

vBoundaries(Layers - 1,1) = (a(Layers - 1)^2)-(b(Layers - 1)^2)*(Radii(Layers))^2;
vBoundaries(Layers - 1,2) = 0.34042; %Air

rhoBoundaries = [ 12.76, 12.17; 9.90, 5.57; 1.02, 1.25*(10^(-3))].*(10^12);

%GENERATING THE RANDOM POINTS ON A SPHERE

M = 12;
N = 10;

%GENERATING RANDOM POINTS ON A SPHERE

Th = rand(M,N).*(2*pi);
Ph = asin((rand(M,N).*2.-1));
I = transpose(1:M);
J = 1 : N;
RandSph = zeros(M,N,3);
RandSph(:,:,1) = Th(:,:);
RandSph(:,:,2) = Ph(:,:);
RandSph(I,J,3) = 1;
[RandX,RandY,RandZ] = sph2cart(RandSph(I,J,1),RandSph(I,J,2),RandSph(I,J,3));

%PREPARING THE INITIAL CONDITIONS FOR THE RAYS

L = 2*((((Radii(Layers))^2)-((D)^2))^0.5);
l = repmat((I-((M+1)/2)).*(L/(M+1)),[1,N]);
x0 = repmat(X0(L,(M+1),I,D), [1,N]);
AlphaRay = zeros(M,N);

Alpha = atan2(D,l);
RandX = cos(Alpha).*(RandX + l) + sin(Alpha).*(RandY + D);
RandY = -sin(Alpha).*(RandX + l) + cos(Alpha).*(RandY + D);

Beta = atan2(RandZ,RandY);
RandY = cos(Beta).*RandY + sin(Beta).*RandZ;
RandZ = -sin(Beta).*RandY + cos(Beta).*RandZ;

RayArrayx0 = zeros(M,N,2^(ACO));
RayArrayAmp = zeros(M,N,2^(ACO));
RayArrayAngDisp = zeros(M,N,2^(ACO));
RayArrayAlpha = zeros(M,N,2^(ACO));

%ASSIGNMENT OF INITIAL CONDITIONS

RayArrayAlpha(:,:,1) = atan2(RandY,RandX - x0);
RayArrayx0(:,:,1) = x0(:,:);
RayArrayx0(:,:,2:2^(ACO)) = NaN;
RayArrayAmp(:,:,1) = 1;
RayArrayPlaneAngle = repmat(Alpha,[1,N,2^(ACO)]);
RayArrayRayAngle = repmat(Beta,[1,1,2^(ACO)]);

IntensityArrayAngDisp = zeros(M,N,2^(ACO));
IntensityArrayAmp = zeros(M,N,2^(ACO));

SOutward100 = zeros(M,N,2^(ACO));
SOutward000 = zeros(M,N,2^(ACO));
SOutward001 = zeros(M,N,2^(ACO));
SInward100 = zeros(M,N,2^(ACO));
SInward000 = zeros(M,N,2^(ACO));
SInward001 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceOutward100 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceOutward000 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceOutward001 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceInward100 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceInward000 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceInward001 = zeros(M,N,2^(ACO));
TIROutward100 = zeros(M,N,2^(ACO));
TIROutward000 = zeros(M,N,2^(ACO));
TIROutward001 = zeros(M,N,2^(ACO));
TIRInward100 = zeros(M,N,2^(ACO));
TIRInward000 = zeros(M,N,2^(ACO));
TIRInward001 = zeros(M,N,2^(ACO));
ThetaOutward100 = zeros(M,N,2^(ACO));
ThetaOutward000 = zeros(M,N,2^(ACO));
ThetaOutward001 = zeros(M,N,2^(ACO));
ThetaInward100 = zeros(M,N,2^(ACO));
ThetaInward000 = zeros(M,N,2^(ACO));
ThetaInward001 = zeros(M,N,2^(ACO));


%PROPAGATION

LeftRays = zeros(M,N,2^(ACO),1);
RightRays = zeros(M,N,2^(ACO),1);
LayerRays = zeros(M,N,2^(ACO),Layers - 1);

LIdOutward = zeros(M,N,2^(ACO));
LIdInward = zeros(M,N,2^(ACO));

LIdInwardHit = zeros(M,N,2^(ACO));
LIdInwardMiss = zeros(M,N,2^(ACO));

for i = 1 : ACO
    for j = 1 : Layers - 1
        SOutward100(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),1,0,0);
        SOutward000(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,0);
        SOutward001(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,1);
        sinAngleOfIncidenceOutward100(:,:,1 : 2^(i - 1)) = GOD(SOutward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);
        sinAngleOfIncidenceOutward000(:,:,1 : 2^(i - 1)) = GOD(SOutward000(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);
        sinAngleOfIncidenceOutward001(:,:,1 : 2^(i - 1)) = GOD(SOutward001(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);
        ThetaOutward100(:,:,1 : 2^(i - 1)) = theta(SOutward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);
        ThetaOutward000(:,:,1 : 2^(i - 1)) = theta(SOutward000(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);
        ThetaOutward001(:,:,1 : 2^(i - 1)) = theta(SOutward001(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);
        SInward100(:,:,1 : 2^(i - 1)) = SR(Radii(j),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),1,0,0);
        SInward000(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,0);
        SInward001(:,:,1 : 2^(i - 1)) = SR(Radii(j),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,1);
        sinAngleOfIncidenceInward100(:,:,1 : 2^(i - 1)) = GOD(SInward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);
        sinAngleOfIncidenceInward000(:,:,1 : 2^(i - 1)) = GOD(SInward000(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);
        sinAngleOfIncidenceInward001(:,:,1 : 2^(i - 1)) = GOD(SInward001(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);
        %sinAngleOfIncidence =                          GOD(SInward001(LIdInwardHit),a(j),b(j), pi - RayArrayAlpha(LIdInwardHit),RayArrayx0(LIdInwardHit),0,0,1);
        %sinAngleOfIncidence =                          GOD(SInwardMiss(LIdInwardMiss),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardMiss),RayArrayx0(LIdInwardMiss),0,0,0);
        ThetaInward100(:,:,1 : 2^(i - 1)) = theta(SInward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);
        ThetaInward000(:,:,1 : 2^(i - 1)) = theta(SInward001(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);
        ThetaInward001(:,:,1 : 2^(i - 1)) = theta(SInward001(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);
        %                                   theta(SInward001(LIdInwardHit),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardHit),RayArrayx0(LIdInwardHit),0,0,1);
        %                                   theta(SInwardMiss(LIdInwardMiss),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardMiss),RayArrayx0(LIdInwardMiss),0,0,0);
        TIROutward100(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward100(:,:,1 : 2^(i - 1))) > 1.000000000000001;
        TIROutward000(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(:,:,1 : 2^(i - 1))) > 1.000000000000001;
        TIROutward001(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward001(:,:,1 : 2^(i - 1)))> 1.000000000000001;
        TIRInward100(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward100(:,:,1 : 2^(i - 1))) > 1.000000000000001;
        TIRInward000(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward000(:,:,1 : 2^(i - 1))) > 1.000000000000001;%miss
        TIRInward001(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward001(:,:,1 : 2^(i - 1))) > 1.000000000000001;%hit
        
        if i == 1
            LIdOutward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) < pi / 2 ) & ( Radii(j) <= RayArrayx0(:,:,1 : 2^(i - 1)) & RayArrayx0(:,:,1 : 2^(i - 1)) < Radii(j + 1) ) & (~LIdOutward(:,:,1 : 2^(i-1)));
            LIdInward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) >= pi / 2 ) & ( Radii(j) <= RayArrayx0(:,:,1 : 2^(i - 1)) & RayArrayx0(:,:,1 : 2^(i - 1)) < Radii(j + 1) ) & (~LIdInward(:,:,1 : 2^(i-1)));
        else %Every ray is on a boundary
            if j == Layers - 1
                LIdOutward(:,:,1 : 2^(i - 1)) = false(M,N,2^(i - 1)); %Outward bound rays cannot exist on the outermost boundary
            else
                LIdOutward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) < pi / 2 ) & ( Radii(j + 1) == RayArrayx0(:,:,1 : 2^(i - 1)) ) & (~LIdOutward(:,:,1 : 2^(i-1)));
            end
            LIdInward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) >= pi / 2 ) & ( Radii(j + 1) == RayArrayx0(:,:,1 : 2^(i - 1)) ) & (~LIdInward(:,:,1 : 2^(i-1)));
        end
        
        %OUTWARD BOUND RAYS
        
        TempId = cat(3,LIdOutward(:,:,1 : 2^(i-1)) & ~TIROutward100(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
        TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
        TempIdTIR = cat(3,LIdOutward(:,:,1 : 2^(i-1)) & TIROutward100(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
        
        %REFRACTED OUTWARD BOUND RAYS
        
        RayArrayAngDisp(TempId) = RayArrayAngDisp(TempId) + ThetaOutward100(TempId);
        RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceOutward100(TempId));
        RayArrayx0(TempId) = Radii(j + 1);
        RayArrayAmp(TempId) = TransmissionCoefficient(j,sinAngleOfIncidenceOutward100(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        
        %REFLECTED OUTWARD BOUND RAYS
        
        RayArrayAlpha(TempIdPrime) = pi - asin(sinAngleOfIncidenceOutward100(TempId));
        RayArrayx0(TempIdPrime) = Radii(j + 1);
        RayArrayAngDisp(TempIdPrime) = RayArrayAngDisp(TempId) + ThetaOutward100(TempId);
        RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,sinAngleOfIncidenceOutward100(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
        RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaOutward100(TempIdTIR);
        RayArrayAlpha(TempIdTIR) = pi - asin(sinAngleOfIncidenceOutward100(TempIdTIR));
        RayArrayx0(TempIdTIR) = Radii(j + 1);
        
        %INWARD BOUND RAYS
        
        %DO THE RAYS GO DOWN A LAYER OR UP A LAYER?
        
        %INWARD HIT CHECK

        LIdInwardHit = (abs(imag(SInward001(:,:,1 : 2^(i - 1)))) < 0.000000000000001) & LIdInward;
        LIdInwardMiss = (~LIdInwardHit) & LIdInward;
        
        TempId = cat(3,LIdInwardMiss(:,:,1 : 2^(i-1)) & ~TIRInward000(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
        TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
        TempIdTIR = cat(3,LIdInwardMiss(:,:,1 : 2^(i-1)) & TIRInward000(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
            
        %RAYS THAT MISS THE INNER LAYER
        
        %REFRACTED OUTWARD BOUND RAYS

        RayArrayAngDisp(TempId) = RayArrayAngDisp(TempId) + ThetaInward000(TempId);
        RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempId) = Radii(j + 1);
        RayArrayAmp(TempId) = TransmissionCoefficient(j,sinAngleOfIncidenceInward000(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));

        %REFLECTED OUTWARD BOUND RAYS

        RayArrayAlpha(TempIdPrime) = pi - asin(sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempIdPrime) = Radii(j + 1);
        RayArrayAngDisp(TempIdPrime) = RayArrayAngDisp(TempId) + ThetaInward000(TempId);
        RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,sinAngleOfIncidenceInward000(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
        RayArrayAngDisp(TempIdTIR) = pi - asin(sinAngleOfIncidenceInward000(TempId));
        RayArrayAlpha(TempIdTIR) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempIdTIR) = Radii(j + 1);

        %RAYS THAT HIT THE INNER LAYER

            %REFRACTED INWARD BOUND RAYS

            if j > 1
                RayArrayAlpha(LIdInwardHit) = asin(((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward001(LIdInwardHit)));
                RayArrayx0(LIdInwardHit) = Radii(j);
                RayArrayAngDisp(LIdInwardHit) = RayArrayAngDisp(LIdInwardHit) + ThetaInward001(LIdInwardHit);
                RayArrayAmp(LIdInwardHit) = TransmissionCoefficient(j,sinAngleOfIncidenceInward001(LIdInwardHit),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));

            %REFLECTED INWARD BOUND RAYS

                LIdInwardHitPrime = circshift(LIdInwardHit,[0,0, 2^(i - 1)]);
                RayArrayAlpha(LIdInwardHitPrime) = 2*pi - asin(sinAngleOfIncidenceInward001(LIdInwardHit));
                RayArrayx0(LIdInwardHitPrime) = Radii(j);
                RayArrayAngDisp(LIdInwardHitPrime) = RayArrayAngDisp(LIdInwardHit) + ThetaInward001(LIdInwardHit);
                RayArrayAmp(LIdInwardHitPrime) = ReflectionCoefficient(j,sinAngleOfIncidenceInward001(LIdInwardHit),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
                
            %TOTALLY REFLECTED INWARD BOUND RAYS
            end
        
        
        
    end
end