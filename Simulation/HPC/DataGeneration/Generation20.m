clear

format long

%INPUTS
D = 20*85*10^3; %THE DISTANCE OF CLOSEST APPROACH TO THE CENTER m
ACO = 5; %IT TAKES 16 ITERATIONS FOR A RAY TO PASS ENTIRELY THROUGH THE MOON FOR UNSIMPLIFIED ATTENUATION. IT TAKES ONLY 8 FOR SIMPLIFIED ATTENUATION.
M = 400; %NUMBER OF POINTS ON THE LINE
N = 12 * 2^8; %NUMBER OF RAYS FROM EACH POINT

%VPREMOON MODEL WITH SIMPLIFIED ATTENUATION

%%{
a = ([8.26955,5.5,3.2,1.0].*10^3).^(1/2); %THESE ARE THE CONSTANT TERMS IN THE PARABOLIC VELOCITY FIELD STRUCTURE, THE FIRST ENTRY CORRESPONDING TO THE CORE km/s
b = ([2.47697*10^(-7),10^(-11),10^(-11),10^(-11)].*10^(-3)).^(1/2); %THESE ARE THE COEFFICIENTS OF THE QUADRATIC TERM IN THE VELOCITY FIELD 1/(s*km)

as = ([4.66221,3.3,1.8,0.5].*10^3).^(1/2); %SHEAR WAVE VELOCITY CONSTANTS
bs = ([1.08519*10^(-7),10^(-11),10^(-11),10^(-11)].*10^(-3)).^(1/2); %NOTE THAT b IS ESSENTIALLY SET TO 0 FOR 3 OUTER LAYERS BECAUSE THEY ARE SO THIN

Radii = [0,1709.1,1725.1,1736.1,1737.1].*10^3; %NOTE THAT WE MUST TAKE 0 AS A LAYER FOR TECHNICAL PURPOSES km

Q = [6750,6750,6750,6750].*10^3; %THE QUALITY FACTOR IN THE DIFFERENT REGIONS

AverageV = [8.07936,5.5,3.2,1.0].*10^3; %THE AVERAGE VELOCITY IN EACH OF THESE REGIONS km/s (INTEGRAL r^2 dr v(r) / INTEGRAL r^2 dr), USED TO CALCULATE ATTENUATION, NOT RAY PROPAGATION

rhoBoundaries = [   2.762       ,   2.762               ;...
                    2.762       ,   2.762               ;...
                    2.762       ,   2.600               ;...
                    2.600       ,   0.000               ].*(10^3);
                    %THE DENSITY OF THE MEDIA AT THEIR DISCONTINUOUS
                    %BOUNDARIES - THE LEFT ENTRY IN EACH ROW CORRESPONDS TO
                    %THE INNER DENSITY AND THE RIGHT CORRESPONDS TO THE
                    %OUTER DENSITY kg/m^3
%}

%VPREMOON MODEL (PLOTTED DATA POINTS AND FIT)

%{
a = ([8.26955,8.26955,8.26955,8.26955,8.26955,5.5,3.2,1.0].*10^3).^(1/2); %THESE ARE THE CONSTANT TERMS IN THE PARABOLIC VELOCITY FIELD STRUCTURE, THE FIRST ENTRY CORRESPONDING TO THE CORE km/s
b = ([2.47697*10^(-7),2.47697*10^(-7),2.47697*10^(-7),2.47697*10^(-7),2.47697*10^(-7),10^(-11),10^(-11),10^(-11)].*10^(-3)).^(1/2); %THESE ARE THE COEFFICIENTS OF THE QUADRATIC TERM IN THE VELOCITY FIELD 1/(s*km)

as = ([4.66221,4.66221,4.66221,4.66221,4.66221,3.3,1.8,0.5].*10^3).^(1/2); %SHEAR WAVE VELOCITY CONSTANTS
bs = ([1.08519*10^(-7),1.08519*10^(-7),1.08519*10^(-7),1.08519*10^(-7),1.08519*10^(-7),10^(-11),10^(-11),10^(-11)].*10^(-3)).^(1/2); %NOTE THAT b IS ESSENTIALLY SET TO 0 FOR 3 OUTER LAYERS BECAUSE THEY ARE SO THIN

Radii = [0,961.7,1231.7,1461.7,1671.7,1709.1,1725.1,1736.1,1737.1].*10^3; %NOTE THAT WE MUST TAKE 0 AS A LAYER FOR TECHNICAL PURPOSES km

Q = [675,1125,3375,9000,6750,6750,6750,6750].*10^3; %THE QUALITY FACTOR IN THE DIFFERENT REGIONS

AverageV = [8.15502,7.96714,7.81709,7.65888,7.56174,5.5,3.2,1.0].*10^3; %THE AVERAGE VELOCITY IN EACH OF THESE REGIONS km/s (INTEGRAL r^2 dr v(r) / INTEGRAL r^2 dr), USED TO CALCULATE ATTENUATION, NOT RAY PROPAGATION

rhoBoundaries = [   3.406       ,   3.406               ;...
                    3.379       ,   3.379               ;...
                    3.350       ,   3.350               ;...
                    3.318       ,   3.318               ;...
                    2.762       ,   2.762               ;...
                    2.762       ,   2.762               ;...
                    2.762       ,   2.600               ;...
                    2.600       ,   0.000               ].*(10^3);
                    %THE DENSITY OF THE MEDIA AT THEIR DISCONTINUOUS
                    %BOUNDARIES - THE LEFT ENTRY IN EACH ROW CORRESPONDS TO
                    %THE INNER DENSITY AND THE RIGHT CORRESPONDS TO THE
                    %OUTER DENSITY kg/m^3
%}

%DENSITY FIT

c = 3.446 * 10^3; %THE CONSTANTS WHICH DESCRIBE THE DENSITY FIELD c - d(r - e)^2 kg/m^3
d = 5*(10^(-8)) * 10^(-3); %THESE ARE USED TO GET THE INITIAL ENERGY OF THE RAYS d=kg/m^5
e = 70.*10^3; %NOTE THE DIFFERENCE IN UNITS BETWEEN c, d, AND e! m

%END INPUTS

Layers = length(Radii);

vBoundaries = zeros(Layers-2,2); %v(:,1) corresponds to the inner velocity on the :th layer, and v(:,2) corresponds to the outer velocity on the :th layer
ReflectionCoefficients = zeros(Layers - 1,2); %NOT SURE WHAT THIS MEANS: these correspond to 4 functions of the angle of incidence
TransmissionCoefficients = zeros(Layers - 1,2); %these correspond to 4 functions of the angle of incidence
for i = 1 : Layers - 2
    vBoundaries(i,1) = (a(i)^2)-(b(i)^2)*(Radii(i+1))^2; % vB(:,1) CORRESPONDS TO THE INNER VELOCITY ON THE :th LAYER, and vB(:,2) CORRESPONDS TO THE OUTER VELOCITY ON THE :th LAYER
    vBoundaries(i,2) = (a(i+1)^2)-(b(i+1)^2)*(Radii(i+1))^2; %outer
end

vBoundaries(Layers - 1,1) = (a(Layers - 1)^2)-(b(Layers - 1)^2)*(Radii(Layers))^2;
vBoundaries(Layers - 1,2) = 0; %AIR/SPACE

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

L = 2*((((Radii(Layers))^2)-((D)^2))^0.5); %LENGTH OF MINIMUM IMPACT PARAMETER MACRO CORD
l = repmat((I-((M+1)/2)).*(L/(M+1)),[1,N]); % PARTITION OF L, THE M+1 PREVENTS THE LAST POINT FROM BEING ON THE SURFACE OF THE MOON. THERE IS NO 0th POINT
x0 = repmat(X0(L,(M+1),I,D), [1,N]); %DISTANCE OF EACH POINT FROM THE ORIGIN (EXTERNALLY DEFINED FUNCTION)
AlphaRay = zeros(M,N); %ANGLE FROM THE CENTER OF THE MOON TO EACH POINT ALONG THE MACRO TRAJECTORY

%WE NOW ROTATE OUR FRAME SUCH THAT THE TAIL OF EACH VECTOR (DEFINING THE INITIAL TRAJECTORY LOCATION) IS IN THE XY PLANE

Alpha = atan2(D,l);
RandXrote = cos(Alpha).*(RandX + l) + sin(Alpha).*(RandY + D);
RandYrote = -sin(Alpha).*(RandX + l) + cos(Alpha).*(RandY + D);
RandZrote = RandZ;

%WE NOW ROTATE OUR FRAME SUCH THAT THE TIP OF EACH VECTOR (DEFINING THE INITIAL TRAJECTORY DIRECTION) IS IN THE XY PLANE

Beta = atan2(RandZrote,RandYrote);
RandX = RandXrote;
RandY = cos(Beta).*RandYrote + sin(Beta).*RandZrote;
RandZ = -sin(Beta).*RandYrote + cos(Beta).*RandZrote;

SourceVelocity = zeros(size(x0));
mu = zeros(size(x0)); %LAMÈ COEFFICIENTS
lambda = zeros(size(x0));
BulkMod = zeros(size(x0));
kappa = zeros(size(x0));

for i = 1 : Layers - 1
    index = (x0 > Radii(i) & x0< Radii(i+1));
    SourceVelocity(index) = a(i).^2 - (b(i).^2).*(x0(index).^2);
    mu(index) = ((as(i).^2 - (bs(i).^2).*(x0(index).^2)).^2).*(c - d.*((x0(index) - e).^2));
    lambda(index) = (((a(i).^2 - (b(i).^2).*(x0(index).^2)).^2) - 2.*((as(i)^2 - (bs(i)^2).*(x0(index).^2)).^2)).*(c - d.*((x0(index) - e).^2));
    BulkMod(index) = lambda(index) + (2/3).*mu(index);
    kappa(index) = (BulkMod(index).^2)./(lambda(index) + 2 .* mu(index));
end

RayArrayx0 = zeros(M,N,2^(ACO));
RayArrayAmp = zeros(M,N,2^(ACO));
RayArraySup = zeros(M,N,2^(ACO));
RayArrayAngDisp = zeros(M,N,2^(ACO));
RayArrayAlpha = zeros(M,N,2^(ACO));
RayArrayTime = zeros(M,N,2^(ACO));

%ASSIGNMENT OF INITIAL CONDITIONS

RayArrayAlpha(:,:,1) = atan2(RandY,RandX - x0);
RayArrayx0(:,:,1) = x0(:,:);
RayArrayx0(:,:,2:2^(ACO)) = NaN;
RayArrayAmp(:,:,1) = 1;
RayArrayPlaneAngle = repmat(Alpha,[1,1,2^(ACO + 1)]);
RayArrayRayAngle = repmat(Beta,[1,1,2^(ACO + 1)]);
kappa = repmat(kappa,[1,1,2^(ACO + 1)]);
SourceVelocity = repmat(SourceVelocity,[1,1,2^(ACO + 1)]);

Density = c - d.*((RayArrayx0(:,:,1) - e).^2);

RayArrayDensity = repmat(Density,[1,1,2^(ACO + 1)]);

OutArrayAngDisp = zeros(M,N,2^(ACO + 1));
OutArrayAmp = zeros(M,N,2^(ACO + 1));
OutArraySup = zeros(M,N,2^(ACO + 1));
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
        else %EVERY RAY IS ON A BOUNDARY
            if j == Layers - 1
                LIdOutward(:,:,1 : 2^(i - 1)) = false(M,N,2^(i - 1)); %OUTWARD BOUND RAYS CANNOT EXIST ON THE OUTERMOST BOUNDARY
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
            TempSup = RayArraySup(TempId);
        
            RayArrayAngDisp(TempId) = TempDisp + ThetaOutward100(TempId);
            RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceOutward100(TempId));
            if i == 1
                RayArrayx0(TempId) = Radii(j + 1);
            elseif j < Layers - 1
                RayArrayx0(TempId) = Radii(j + 2);
            end
            RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceOutward100(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
            RayArraySup(TempId) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
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
            RayArraySup(TempIdPrime) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
            RayArrayTime(TempIdPrime) = TempTime + SOutward100(TempId) - SXOutward(TempId);
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
            RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaOutward100(TempIdTIR);
            RayArrayAlpha(TempIdTIR) = pi - asin(sinAngleOfIncidenceOutward100(TempIdTIR));
            RayArrayx0(TempIdTIR) = Radii(j + 1);
            RayArraySup(TempIdTIR) = QualityFactor(j,SOutward100(TempIdTIR) - SXOutward(TempIdTIR),AverageV,Q,RayArraySup(TempIdTIR));
            RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SOutward100(TempIdTIR) - SXOutward(TempIdTIR);
        
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
        TempSup = RayArraySup(TempId);
        
        RayArrayAngDisp(TempId) = TempDisp + ThetaInward000(TempId);
        RayArrayAlpha(TempId) = asin((vBoundaries(j,2)/(vBoundaries(j,1))).*sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempId) = Radii(j + 1);
        RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceInward000(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        RayArraySup(TempId) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
        RayArrayTime(TempId) = TempTime + SInward000(TempId) - SXInward(TempId);
        
        %REFLECTED OUTWARD BOUND RAYS (THESE ARE THE RAYS WHICH CAN BOUNCE ALONG THE SURFACE)
        
        RayArrayAngDisp(TempIdPrime) = TempDisp + ThetaInward000(TempId);
        RayArrayAlpha(TempIdPrime) = pi - asin(sinAngleOfIncidenceInward000(TempId));
        RayArrayx0(TempIdPrime) = Radii(j + 1);
        RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,TempAmp,sinAngleOfIncidenceInward000(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
        RayArraySup(TempIdPrime) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
        RayArrayTime(TempIdPrime) = TempTime + SInward000(TempId) - SXInward(TempId);
        
        %TOTALLY REFLECTED OUTWARD BOUND RAYS
        
        RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaInward000(TempIdTIR);
        RayArrayAlpha(TempIdTIR) = pi - asin(sinAngleOfIncidenceInward000(TempIdTIR));
        RayArrayx0(TempIdTIR) = Radii(j + 1);
        RayArraySup(TempIdTIR) = QualityFactor(j,SOutward100(TempIdTIR) - SXOutward(TempIdTIR),AverageV,Q,RayArraySup(TempIdTIR));
        RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SInward000(TempIdTIR) - SXInward(TempIdTIR);
        
        %RAYS THAT HIT THE INNER LAYER

            if j > 1
                
            %REFRACTED INWARD BOUND RAYS

                TempId = cat(3,LIdInwardHit(:,:,1 : 2^(i-1)) & ~TIRInward001(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
                TempIdPrime = circshift(TempId,[0,0, 2^(i - 1)]);
                TempIdTIR = cat(3,LIdInwardHit(:,:,1 : 2^(i-1)) & TIRInward001(:,:,1 : 2^(i-1)),false(M,N,2^(i-1)));
                
                TempDisp = RayArrayAngDisp(TempId);
                TempTime = RayArrayTime(TempId);
                TempAmp = RayArrayAmp(TempId);
                TempSup = RayArraySup(TempId);
                
                if i == 1
                    RayArrayAlpha(TempId) = pi - asin(((vBoundaries(j - 1,1)/(vBoundaries(j - 1,2))).*sinAngleOfIncidenceInward001(TempId)));
                else
                    RayArrayAlpha(TempId) = pi - asin(((vBoundaries(j,1)/(vBoundaries(j,2))).*sinAngleOfIncidenceInward001(TempId)));
                end
                RayArrayx0(TempId) = Radii(j);
                RayArrayAngDisp(TempId) = TempDisp + ThetaInward001(TempId);
                RayArrayAmp(TempId) = TransmissionCoefficient(j,TempAmp,sinAngleOfIncidenceInward001(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
                RayArraySup(TempId) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
                RayArrayTime(TempId) = TempTime + SInward001(TempId) - SXInward(TempId);
                
            %REFLECTED INWARD BOUND RAYS
                
                RayArrayAlpha(TempIdPrime) = asin(sinAngleOfIncidenceInward001(TempId));
                RayArrayx0(TempIdPrime) = Radii(j);
                RayArrayAngDisp(TempIdPrime) = TempDisp + ThetaInward001(TempId);
                RayArrayAmp(TempIdPrime) = ReflectionCoefficient(j,TempAmp,sinAngleOfIncidenceInward001(TempId),RayArrayAlpha(TempId),vBoundaries(j,1),vBoundaries(j,2),rhoBoundaries(j,1),rhoBoundaries(j,2));
                RayArraySup(TempIdPrime) = QualityFactor(j,SOutward100(TempId) - SXOutward(TempId),AverageV,Q,TempSup);
                RayArrayTime(TempIdPrime) = TempTime + SInward001(TempId) - SXInward(TempId);
                
            %TOTALLY REFLECTED INWARD BOUND RAYS
            
                RayArrayx0(TempIdTIR) = Radii(j);
                RayArrayAlpha(TempIdTIR) = asin(sinAngleOfIncidenceInward001(TempIdTIR));
                RayArrayAngDisp(TempIdTIR) = RayArrayAngDisp(TempIdTIR) + ThetaInward001(TempIdTIR);
                RayArraySup(TempIdTIR) = QualityFactor(j,SOutward100(TempIdTIR) - SXOutward(TempIdTIR),AverageV,Q,RayArraySup(TempIdTIR));
                RayArrayTime(TempIdTIR) = RayArrayTime(TempIdTIR) + SInward001(TempIdTIR) - SXInward(TempIdTIR);
                
            end
    end
    
    IdSurface = RayArrayx0 == Radii(Layers);
    IdSurfacePrime = cat(3,false(M,N,2^(i - 1)),IdSurface);

    OutArrayAngDisp(IdSurfacePrime) = RayArrayAngDisp(IdSurface);
    OutArrayAmp(IdSurfacePrime) = RayArrayAmp(IdSurface);
    OutArraySup(IdSurfacePrime) = RayArraySup(IdSurface);
    OutArrayTime(IdSurfacePrime) = RayArrayTime(IdSurface);
    
end
REF = OutArrayAmp(OutArrayTime > 0 & OutArrayAmp > 0);
ATT= OutArraySup(OutArrayTime > 0 & OutArrayAmp > 0);
OutTime = OutArrayTime(OutArrayTime > 0 & OutArrayAmp > 0);
Density = RayArrayDensity(OutArrayTime > 0 & OutArrayAmp > 0);
kappa = kappa(OutArrayTime > 0 & OutArrayAmp> 0 );
velocity = SourceVelocity(OutArrayTime > 0 & OutArrayAmp > 0);

OutAng = OutArrayAngDisp(OutArrayTime > 0 & OutArrayAmp > 0);
OutRay = RayArrayRayAngle(OutArrayTime > 0 & OutArrayAmp > 0);
OutPla = RayArrayPlaneAngle(OutArrayTime > 0 & OutArrayAmp > 0);

[X,Y,Z] = sph2cart(OutAng,0,1);
Xrote = X;
Yrote = cos(OutRay).*Y - sin(OutRay).*Z;
Zrote = sin(OutRay).*Y + cos(OutRay).*Z;
Xfinal = cos(OutPla).*Xrote - sin(OutPla).*Yrote;
Yfinal = sin(OutPla).*Xrote + cos(OutPla).*Yrote;
Zfinal = Zrote;

outputmat = cat(2,Xfinal,Yfinal,Zfinal,OutTime,Density,REF,ATT,kappa,velocity);
outputnum = [M,N,L,D];

save Mat20.mat outputmat;
save Num20.mat outputnum;