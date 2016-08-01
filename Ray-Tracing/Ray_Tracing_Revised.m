clear

format long



%INPUTS

M = 1; %NUMBER OF POINTS ON THE LINE
N = 5000; %NUMBER OF RAYS FROM EACH POINT

a = [11.2633,10.6835,15.4938].^(1/2); %THESE ARE THE CONSTANT TERMS IN THE PARABOLIC VELOCITY FIELD STRUCTURE, THE FIRST ENTRY CORRESPONDING TO THE CORE
b = [1.55946*10^(-7),2.11853*10^(-7),1.40336*10^(-7)].^(1/2); %THESE ARE THE COEFFICIENTS OF THE QUADRATIC TERM IN THE VELOCITY FIELD

D = 6370; %THE DISTANCE OF CLOSEST APPROACH TO THE CENTER
Radii = [0,1221,3480,6371]; %NOTE THAT WE MUST TAKE 0 AS A LAYER FOR TECHNICAL PURPOSES
ACO = 8; %IT TAKES 5 ITERATIONS FOR A RAY TO PASS ENTIRELY THROUGH THE EARTH

%Q = ; %THE QUALITY FACTOR AS A FUNCTION OF FREQUENCY, TIME, AND SO ON... TBD
%A = ; %THE AMPLITUDE OF THE WAVE AS A FUNCTION OF FREQUENCY

rhoBoundaries = [ 12.76, 12.17; 9.90, 5.57; 1.02, 1.25*(10^(-3))].*(10^12); %THE DENSITY OF THE MEDIA AT THEIR DISCONTINUOUS BOUNDARIES - THE LEFT ENTRY IN EACH ROW CORRESPONDS TO THE INNER DENSITY AND THE RIGHT CORRESPONDS TO THE OUTER DENSITY

%END INPUTS









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
RandXrote = cos(Alpha).*(RandX + l) + sin(Alpha).*(RandY + D);
RandYrote = -sin(Alpha).*(RandX + l) + cos(Alpha).*(RandY + D);
RandZrote = RandZ;

Beta = atan2(RandZrote,RandYrote);
RandX = RandXrote;
RandY = cos(Beta).*RandYrote + sin(Beta).*RandZrote;
RandZ = -sin(Beta).*RandYrote + cos(Beta).*RandZrote;

RayArrayx0 = zeros(M,N,2^(ACO));
RayArrayAmp = zeros(M,N,2^(ACO));
RayArrayAngDisp = zeros(M,N,2^(ACO));
RayArrayAlpha = zeros(M,N,2^(ACO));
RayArrayTime = zeros(M,N,2^(ACO));

%ASSIGNMENT OF INITIAL CONDITIONS
%%{
RayArrayAlpha(:,:,1) = atan2(RandY,RandX - x0);
RayArrayx0(:,:,1) = x0(:,:);
RayArrayx0(:,:,2:2^(ACO)) = NaN;
RayArrayAmp(:,:,1) = 1;
RayArrayPlaneAngle = repmat(Alpha,[1,1,2^(ACO + 1)]);
RayArrayRayAngle = repmat(Beta,[1,1,2^(ACO + 1)]);
%}
%TESTING
%{
                RayArrayAlpha(:,:,1) = pi;
                RayArrayx0(:,:,1) = 6370;
                RayArrayx0(:,:,2:2^(ACO)) = NaN;
                RayArrayAmp(:,:,1) = 1;
                RayArrayPlaneAngle = zeros(M,N,2^(ACO+1));
                RayArrayRayAngle = zeros(M,N,2^(ACO+1));
%}
%END TESTING

OutArrayAngDisp = zeros(M,N,2^(ACO + 1));
OutArrayAmp = zeros(M,N,2^(ACO + 1));
OutArrayTime = zeros(M,N,2^(ACO + 1));

SOutward100 = zeros(M,N,2^(ACO));
SInward000 = zeros(M,N,2^(ACO));
SInward001 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceOutward100 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceInward000 = zeros(M,N,2^(ACO));
sinAngleOfIncidenceInward001 = zeros(M,N,2^(ACO));
TIROutward100 = false(M,N,2^(ACO));
TIRInward000 = false(M,N,2^(ACO));
TIRInward001 = false(M,N,2^(ACO));
ThetaOutward100 = zeros(M,N,2^(ACO));
ThetaInward000 = zeros(M,N,2^(ACO));
ThetaInward001 = zeros(M,N,2^(ACO));
SXOutward = zeros(M,N,2^(ACO));
SXInward = zeros(M,N,2^(ACO));

%PROPAGATION

LIdOutward = false(M,N,2^(ACO));
LIdInward = false(M,N,2^(ACO));

LIdInwardHit = false(M,N,2^(ACO));
LIdInwardMiss = false(M,N,2^(ACO));

sx = zeros(M,N);


for i = 1 : ACO
    for j = 1 : Layers - 1
        if or (j < Layers - 1 , i == 1)
            if i == 1
                SOutward100(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),1,0,0);%KEEP
            elseif j < Layers - 1
                SOutward100(:,:,1 : 2^(i - 1)) = SR(Radii(j + 2),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),1,0,0);%KEEP
            end
            sinAngleOfIncidenceOutward100(:,:,1 : 2^(i - 1)) = GOD(SOutward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);%KEEP
            ThetaOutward100(:,:,1 : 2^(i - 1)) = theta(SOutward100(:,:,1 : 2^(i - 1)),a(j),b(j),RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),1,0,0);%KEEP
            TIROutward100(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceOutward100(:,:,1 : 2^(i - 1))) > 1.000000000000001;%KEEP
            SXOutward(:,:,1 : 2^(i - 1)) = SX(a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),RayArrayAlpha(:,:,1 : 2^(i - 1)),1,0,0);
        end
        
        SInward000(:,:,1 : 2^(i - 1)) = SR(Radii(j + 1),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,0);%KEEP
        sinAngleOfIncidenceInward000(:,:,1 : 2^(i - 1)) = GOD(SInward000(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);%KEEP
        ThetaInward000(:,:,1 : 2^(i - 1)) = theta(SInward000(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,0);%KEEP
        TIRInward000(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(:,:,1 : 2^(i - 1))) > 1.000000000000001;%KEEP
        SXInward(:,:,1 : 2^(i - 1)) = SX(a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,0);
        
        if j > 1
            SInward001(:,:,1 : 2^(i - 1)) = SR(Radii(j),a(j),b(j),RayArrayx0(:,:,1 : 2^(i - 1)),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),0,0,1);%KEEP
            sinAngleOfIncidenceInward001(:,:,1 : 2^(i - 1)) = GOD(SInward001(:,:,1 : 2^(i - 1)),a(j),b(j), pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);%KEEP
            ThetaInward001(:,:,1 : 2^(i - 1)) = theta(SInward001(:,:,1 : 2^(i - 1)),a(j),b(j),pi - RayArrayAlpha(:,:,1 : 2^(i - 1)),RayArrayx0(:,:,1 : 2^(i - 1)),0,0,1);%KEEP
            TIRInward001(:,:,1 : 2^(i - 1)) = abs((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward001(:,:,1 : 2^(i - 1))) > 1.000000000000001;%KEEP
        end
        
        if i == 1
            LIdOutward(:,:,1 : 2^(i - 1)) = ((( RayArrayAlpha(:,:,1 : 2^(i - 1)) < pi / 2 ) & ( Radii(j) <= RayArrayx0(:,:,1 : 2^(i - 1)) & RayArrayx0(:,:,1 : 2^(i - 1)) < Radii(j + 1) )) & (~LIdOutward(:,:,1 : 2^(i-1)))) & (~LIdInward(:,:,1 : 2^(i - 1)) );
            LIdInward(:,:,1 : 2^(i - 1)) = ((( RayArrayAlpha(:,:,1 : 2^(i - 1)) >= pi / 2 ) & ( Radii(j) <= RayArrayx0(:,:,1 : 2^(i - 1)) & RayArrayx0(:,:,1 : 2^(i - 1)) < Radii(j + 1) )) & (~LIdInward(:,:,1 : 2^(i-1)))) & (~LIdOutward(:,:,1 : 2^(i - 1)) );
        else %Every ray is on a boundary
            if j == Layers - 1
                LIdOutward(:,:,1 : 2^(i - 1)) = (( RayArrayAlpha(:,:,1 : 2^(i - 1)) < pi / 2 ) & ( Radii(j) <= RayArrayx0(:,:,1 : 2^(i - 1)) & RayArrayx0(:,:,1 : 2^(i - 1)) < Radii(j + 1) )) & (~LIdOutward(:,:,1 : 2^(i-1)));
                RayArrayx0(LIdOutward) = NaN;
                LIdOutward(:,:,1 : 2^(i - 1)) = false(M,N,2^(i - 1)); %Outward bound rays cannot exist on the outermost boundary
            else
                LIdOutward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) < pi / 2 ) & ( Radii(j + 1) == RayArrayx0(:,:,1 : 2^(i - 1)) ) & (~LIdOutward(:,:,1 : 2^(i-1)));
            end
            LIdInward(:,:,1 : 2^(i - 1)) = ( RayArrayAlpha(:,:,1 : 2^(i - 1)) >= pi / 2 ) & ( Radii(j + 1) == RayArrayx0(:,:,1 : 2^(i - 1)) ) & (~LIdInward(:,:,1 : 2^(i-1)));
        end
        
        %OUTWARD BOUND RAYS
        
        if or (j < Layers - 1 , i == 1)
        
            TempId = cat(3,LIdOutward(:,:,1 : 2^(i-1)) & ~TIROutward100(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
            TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
            TempIdTIR = cat(3,LIdOutward(:,:,1 : 2^(i-1)) & TIROutward100(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
        
        %REFRACTED OUTWARD BOUND RAYS
        
            TempDisp = RayArrayAngDisp(TempId);
            TempTime = RayArrayTime(TempId);
            TempAmp = RayArrayAmp(TempId);
        
            RayArrayAngDisp(TempId) = TempDisp + ThetaOutward100(TempId);
            RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceOutward100(TempId));
            if i == 1
                RayArrayx0(TempId) = Radii(j + 1);
            elseif j < Layers - 1
                RayArrayx0(TempId) = Radii(j + 2);
            end
            RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceOutward100(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
            RayArrayTime(TempId) = TempTime + SOutward100(TempId) - SXOutward(TempId);
        
        %REFLECTED OUTWARD BOUND RAYS
        
            RayArrayAngDisp(TempIdPrime) = TempDisp + ThetaOutward100(TempId);
            RayArrayAlpha(TempIdPrime) = pi - asin(sinAngleOfIncidenceOutward100(TempId));
            if i == 1
                RayArrayx0(TempIdPrime) = Radii(j + 1);
            elseif j < Layers - 1
                RayArrayx0(TempIdPrime) = Radii(j + 2);
            end
            RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,TempAmp,sinAngleOfIncidenceOutward100(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
            RayArrayTime(TempIdPrime) = TempTime + SOutward100(TempId) - SXOutward(TempId);
            
        %[nnz(SOutward100(TempId) - SXOutward(TempId) < -0.000000000000002),i,j]
        %[nnz(RayArrayTime(TempId) < -0.000000000000002)]
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
            RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaOutward100(TempIdTIR);
            RayArrayAlpha(TempIdTIR) = pi - asin(sinAngleOfIncidenceOutward100(TempIdTIR));
            RayArrayx0(TempIdTIR) = Radii(j + 1);
            RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SOutward100(TempIdTIR) - SXOutward(TempIdTIR);
            
        %[nnz(SOutward100(TempIdTIR) - SXOutward(TempIdTIR) < -0.000000000000002),i,j]
        %[nnz(RayArrayTime(TempIdTIR) < -0.000000000000002)]
        
        end
        
        %INWARD BOUND RAYS
        
        %DO THE RAYS GO DOWN A LAYER OR UP A LAYER?
        
        %INWARD HIT CHECK
        
        if j > 1
            LIdInwardHit = (abs(imag(SInward001(:,:,1 : 2^(i-1)))) < 0.000000000000001) & LIdInward(:,:,1 : 2^(i-1));
            LIdInwardMiss = ((~LIdInwardHit(:,:,1 : 2^(i-1)))) & LIdInward(:,:,1 : 2^(i-1));
        else
            LIdInwardHit = false(M,N, 2^(i - 1));
            LIdInwardMiss = LIdInward;
        end
        
        TempId = cat(3,(LIdInwardMiss(:,:,1 : 2^(i-1)) & (~TIRInward000(:,:,1 : 2^(i-1)))),false(M,N,2^(i-1)));
        TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
        TempIdTIR = cat(3,(LIdInwardMiss(:,:,1 : 2^(i-1)) & TIRInward000(:,:,1 : 2^(i-1))),false(M,N,2^(i-1)));
            
        %RAYS THAT MISS THE INNER LAYER
        
        %REFRACTED OUTWARD BOUND RAYS
        
        TempDisp = RayArrayAngDisp(TempId);
        TempTime = RayArrayTime(TempId);
        TempAmp = RayArrayAmp(TempId);
        
        RayArrayAngDisp(TempId) = TempDisp + ThetaInward000(TempId);
        RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempId) = Radii(j + 1);
        RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceInward000(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        RayArrayTime(TempId) = TempTime + SInward000(TempId) - SXInward(TempId);
        
        %REFLECTED OUTWARD BOUND RAYS (These are the rays which can bounce
        %along the surface)
        
        RayArrayAngDisp(TempIdPrime) = TempDisp + ThetaInward000(TempId);
        RayArrayAlpha(TempIdPrime) = pi - asin(sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempIdPrime) = Radii(j + 1);
        RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,TempAmp,sinAngleOfIncidenceInward000(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        RayArrayTime(TempIdPrime) = TempTime + SInward000(TempId) - SXInward(TempId);
        
        %[nnz(SInward000(TempId) - SXInward(TempId) < -0.000000000000002),i,j]
        %[nnz(RayArrayTime(TempId) < -0.000000000000002)]
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
        RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaInward000(TempIdTIR);
        RayArrayAlpha(TempIdTIR) = pi - asin(sinAngleOfIncidenceInward000(TempIdTIR));
        RayArrayx0(TempIdTIR) = Radii(j + 1);
        RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SInward000(TempIdTIR) - SXInward(TempIdTIR);
        
        %[nnz(SInward000(TempIdTIR) - SXInward(TempIdTIR) < -0.000000000000002),i,j]
        %[nnz(RayArrayTime(TempIdTIR) < -0.000000000000002)]
        
        %RAYS THAT HIT THE INNER LAYER

            if j > 1
                
            %REFRACTED INWARD BOUND RAYS

                TempId = cat(3,LIdInwardHit(:,:,1 : 2^(i-1)) & ~TIRInward001(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
                TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
                TempIdTIR = cat(3,LIdInwardHit(:,:,1 : 2^(i-1)) & TIRInward001(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
                
                TempDisp = RayArrayAngDisp(TempId);
                TempTime = RayArrayTime(TempId);
                TempAmp = RayArrayAmp(TempId);
                
                RayArrayAlpha(TempId) = pi - asin(((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward001(TempId)));
                RayArrayx0(TempId) = Radii(j);
                RayArrayAngDisp(TempId) = TempDisp + ThetaInward001(TempId);
                RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceInward001(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
                RayArrayTime(TempId) = TempTime + SInward001(TempId) - SXInward(TempId);
                
            %REFLECTED INWARD BOUND RAYS
                
                RayArrayAlpha(TempIdPrime) = asin(sinAngleOfIncidenceInward001(TempId));
                RayArrayx0(TempIdPrime) = Radii(j);
                RayArrayAngDisp(TempIdPrime) = TempDisp + ThetaInward001(TempId);
                RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,TempAmp,sinAngleOfIncidenceInward001(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
                RayArrayTime(TempIdPrime) = TempTime + SInward001(TempId) - SXInward(TempId);
                
                %[nnz(SInward001(TempId) - SXInward(TempId) < -0.000000000000002),i,j]
                %[nnz(RayArrayTime(TempId) < -0.000000000000002)]
                
            %TOTALLY REFLECTED INWARD BOUND RAYS
            
                RayArrayx0(TempIdTIR) = Radii(j);
                RayArrayAlpha(TempIdTIR) = asin(sinAngleOfIncidenceInward001(TempIdTIR));
                RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaInward001(TempIdTIR);
                RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SInward001(TempIdTIR) - SXInward(TempIdTIR);
                
                %[nnz(SInward001(TempIdTIR) - SXInward(TempIdTIR) < -0.000000000000002),i,j]
                %[nnz(RayArrayTime(TempIdTIR) < -0.000000000000002)]
                
            end
    end
    
    IdSurface = RayArrayx0 == Radii(Layers);
    IdSurfacePrime = cat(3,false(M,N,2^(i - 1)),IdSurface);

    OutArrayAngDisp(IdSurfacePrime) = RayArrayAngDisp(IdSurface);
    OutArrayAmp(IdSurfacePrime) = RayArrayAmp(IdSurface);
    OutArrayTime(IdSurfacePrime) = RayArrayTime(IdSurface);
    
end

OutAmp = OutArrayAmp(OutArrayTime > 0);
OutAng = OutArrayAngDisp(OutArrayTime > 0);
OutTim = OutArrayTime(OutArrayTime > 0);
OutRay = RayArrayRayAngle(OutArrayTime > 0);
OutPla = RayArrayPlaneAngle(OutArrayTime > 0);

[X,Y,Z] = sph2cart(OutAng,0,1);
Xrote = X;
Yrote = cos(OutRay).*Y - sin(OutRay).*Z;
Zrote = sin(OutRay).*Y + cos(OutRay).*Z;
Xfinal = cos(OutPla).*Xrote - sin(OutPla).*Yrote;
Yfinal = sin(OutPla).*Xrote + cos(OutPla).*Yrote;
Zfinal = Zrote;
%%{
figure
t0 = 1.5*10^3;
tf = 2*10^3;
scatter3(Xfinal(t0 < OutTim & OutTim < tf),Yfinal(t0 < OutTim & OutTim < tf),Zfinal(t0 < OutTim & OutTim < tf),'.')
%}
%{
figure
scatter3(Xfinal,Yfinal,Zfinal,20,'filled')
%}