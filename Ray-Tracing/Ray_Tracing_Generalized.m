clear

format long

a = [11.2633,10.6835,15.4938].^(1/2); %first entry corresponds to the core
b = [1.55946*10^(-7),2.11853*10^(-7),1.40336*10^(-7)].^(1/2); %first entry corresponds to the core
WeaknessOfNeglectedRays = 10^(-1); %input('What fraction of the initial ray energy must a ray contain to not be thrown out? ');

D = 1500;
Radii = [0,1221,3480,6371];
AbsoluteCutoff = 4;

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

%GENERATING THE RANDOM POINTS ON A SPHERE

M = 1;
N = 100000;
Th = rand(M,N).*(2*pi);
Ph = asin((rand(M,N).*2.-1));
I = 1:M;
J = 1:N;
RandSph = zeros(M,N,3);
RandSph(:,:,1) = Th(:,:);
RandSph(:,:,2) = Ph(:,:);
RandSph(I,J,3) = 1;
[RandX,RandY,RandZ] = sph2cart(RandSph(I,J,1),RandSph(I,J,2),RandSph(I,J,3));

%PREPARING THE INITIAL CONDITIONS FOR THE RAYS

L = 2*((((Radii(Layers))^2)-((D)^2))^0.5);
l = (I-((M+1)/2)).*(L/(M+1));
Alpha = atan2(D,transpose(l));
x0 = transpose(repmat(X0(L,(M+1),I,D), [N,1]));
Beta = zeros(M,N);
AlphaRay = zeros(M,N);
for i = 1 : M
    RotMatrixAlpha = [cos(Alpha(i)),sin(Alpha(i)),0;-sin(Alpha(i)),cos(Alpha(i)),0;0,0,1];
    for j = 1 : N
        RotRand = RotMatrixAlpha*[RandX(i,j)+l(i);RandY(i,j)+D;RandZ(i,j)];
        Beta(i,j) = atan2(RotRand(3),RotRand(2));
        RotMatrixBeta = [1,0,0;0,cos(Beta(i,j)),sin(Beta(i,j));0,-sin(Beta(i,j)),cos(Beta(i,j))];
        Vec = RotMatrixBeta*RotRand;
        AlphaRay(i,j) = atan2(Vec(2),Vec(1) - x0(i,j));
    end
end
RayArrayx0 = zeros(M,N,2^(AbsoluteCutoff));
RayArrayAmp = zeros(M,N,2^(AbsoluteCutoff));
RayArrayAngDisp = zeros(M,N,2^(AbsoluteCutoff));
RayArrayAlpha = zeros(M,N,2^(AbsoluteCutoff));

RayArrayx0(:,:,1) = x0(:,:);
RayArrayx0(:,:,2:2^(AbsoluteCutoff)) = NaN;
RayArrayAlpha(:,:,1) = AlphaRay(:,:);
RayArrayAmp(:,:,1) = 1;
RayArrayPlaneAngle = repmat(Alpha,[1,N,2^(AbsoluteCutoff)]);
RayArrayRayAngle = repmat(Beta,[1,1,2^(AbsoluteCutoff)]);

IntensityArrayAngDisp = zeros(M,N,2^(AbsoluteCutoff));
IntensityArrayAmp = zeros(M,N,2^(AbsoluteCutoff));

%PROPAGATION

LeftRays = zeros(M,N,2^(AbsoluteCutoff),1);
RightRays = zeros(M,N,2^(AbsoluteCutoff),1);
LayerRays = zeros(M,N,2^(AbsoluteCutoff),Layers - 1);

LIdOutward = zeros(M,N,2^(AbsoluteCutoff));
LIdInward = zeros(M,N,2^(AbsoluteCutoff));

for Iter = 1 : AbsoluteCutoff
    for LayerId = 1 : Layers - 1
        if Iter == 1
            LIdOutward = ( RayArrayAlpha(:,:,:) < pi / 2 ) & ( Radii(LayerId) <= RayArrayx0(:,:,:) & RayArrayx0(:,:,:) < Radii(LayerId + 1) ) & (~LIdOutward);
            LIdInward = ( RayArrayAlpha(:,:,:) >= pi / 2 ) & ( Radii(LayerId) <= RayArrayx0(:,:,:) & RayArrayx0(:,:,:) < Radii(LayerId + 1) ) & (~LIdInward);
        else %Every ray is on a boundary
            LIdOutward = ( RayArrayAlpha(:,:,:) < pi / 2 ) & ( Radii(LayerId) == RayArrayx0(:,:,:) ) & (~LIdOutward);
            LIdInward = ( RayArrayAlpha(:,:,:) >= pi / 2 ) & ( Radii(LayerId) == RayArrayx0(:,:,:) ) & (~LIdInward);
        end
        
        %OUTWARD BOUND RAYS
        
            SOutward = SR(Radii(LayerId + 1),a(LayerId),b(LayerId),RayArrayx0(LIdOutward),RayArrayAlpha(LIdOutward),1,0,0);
        
        %REFRACTED OUTWARD BOUND RAYS
        
        sinAngleOfIncidence = GOD(SOutward,a(LayerId),b(LayerId),RayArrayAlpha(LIdOutward),RayArrayx0(LIdOutward),1,0,0);
        Theta = theta(SOutward,a(LayerId),b(LayerId),RayArrayAlpha(LIdOutward),RayArrayx0(LIdOutward),1,0,0);
        if LayerId == Layers - 1
            IntensityArrayAngDisp(LIdOutward) = RayArrayAngDisp(LIdOutward) + Theta;
            IntensityArrayAmp(LIdOutward) = TransmissionCoefficient(sinAngleOfIncidence);
        else
            RayArrayAngDisp(LIdOutward) = RayArrayAngDisp(LIdOutward) + Theta;
            RayArrayAlpha(LIdOutward) = asin((vBoundaries(LayerId,2)/(vBoundaries(LayerId,1))).*sinAngleOfIncidence);
            RayArrayx0(LIdOutward) = Radii(LayerId + 1);
            RayArrayAmp(LIdOutward) = TransmissionCoefficient(sinAngleOfIncidence);
        end

        %REFLECTED OUTWARD BOUND RAYS
        
            LIdOutwardPrime = circshift(LIdOutward,[0,0, 2^(Iter - 1)]);
            RayArrayAlpha(LIdOutwardPrime) = pi - asin(sinAngleOfIncidence);
            RayArrayx0(LIdOutwardPrime) = Radii(LayerId + 1);
            RayArrayAngDisp(LIdOutwardPrime) = RayArrayAngDisp(LIdOutward) + Theta;
            RayArrayAmp(LIdOutwardPrime) = ReflectionCoefficient(sinAngleOfIncidence);
            
        %POTENTIALLY INWARD BOUND RAYS
        
            %INWARD HIT CHECK
            SInwardHit = SR(Radii(LayerId),a(LayerId),b(LayerId),RayArrayx0(:,:,:), pi - RayArrayAlpha(:,:,:),0,0,1);
            SInwardMiss = SR(Radii(LayerId + 1),a(LayerId),b(LayerId),RayArrayx0(:,:,:), pi - RayArrayAlpha(:,:,:),0,0,0);
            LIdInwardHit = (imag(SInwardHit(:,:,:)) == 0) & LIdInward;
            LIdInwardMiss = (~LIdInwardHit) & LIdInward;
            
        %RAYS THAT MISS THE INNER LAYER
        
            %REFRACTED OUTWARD BOUND RAYS
            sinAngleOfIncidence = GOD(SInwardMiss(LIdInwardMiss),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardMiss),RayArrayx0(LIdInwardMiss),0,0,0);
            Thetaprime = theta(SInwardMiss(LIdInwardMiss),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardMiss),RayArrayx0(LIdInwardMiss),0,0,0);
            if LayerId == Layers - 1
                IntensityArrayAngDisp(LIdInwardMiss) = RayArrayAngDisp(LIdInwardMiss) + Thetaprime;
                IntensityArrayAmp(LIdInwardMiss) = TransmissionCoefficient(sinAngleOfIncidence);
            else
                RayArrayAngDisp(LIdInwardMiss) = RayArrayAngDisp(LIdInwardMiss) + Thetaprime;
                RayArrayAlpha(LIdInwardMiss) = asin((vBoundaries(LayerId,2)/(vBoundaries(LayerId,1))).*sinAngleOfIncidence);
                RayArrayx0(LIdInwardMiss) = Radii(LayerId + 1);
                RayArrayAmp(LIdInwardMiss) = TransmissionCoefficient(sinAngleOfIncidence);
            end

            %REFLECTED OUTWARD BOUND RAYS
            
                LIdInwardMissPrime = circshift(LIdInwardMiss,[0,0, 2^(Iter - 1)]);
                RayArrayAlpha(LIdInwardMissPrime) = pi - asin(sinAngleOfIncidence);
                RayArrayx0(LIdInwardMissPrime) = Radii(LayerId + 1);
                RayArrayAngDisp(LIdInwardMissPrime) = RayArrayAngDisp(LIdInwardMissPrime) + Thetaprime;
                RayArrayAmp(LIdInwardMissPrime) = ReflectionCoefficient(sinAngleOfIncidence);
                
        %RAYS THAT HIT THE INNER LAYER
        
                %REFRACTED INWARD BOUND RAYS
                
                if LayerId > 1
                    sinAngleOfIncidence = GOD(SInwardHit(LIdInwardHit),a(LayerId),b(LayerId),RayArrayAlpha(LIdInwardHit),RayArrayx0(LIdInwardHit),0,0,1);
                    Thetaprimeprime = theta(SInwardHit(LIdInwardHit),a(LayerId),b(LayerId), pi - RayArrayAlpha(LIdInwardHit),RayArrayx0(LIdInwardHit),0,0,1);
                    RayArrayAlpha(LIdInwardHit) = asin(((vBoundaries(LayerId,1)/(vBoundaries(LayerId,2))).*sinAngleOfIncidence));
                    RayArrayx0(LIdInwardHit) = Radii(LayerId);
                    RayArrayAngDisp(LIdInwardHit) = RayArrayAngDisp(LIdInwardHit) + Thetaprimeprime;
                    RayArrayAmp(LIdInwardHit) = TransmissionCoefficient(sinAngleOfIncidence);

                %REFLECTED INWARD BOUND RAYS
                
                    LIdInwardHitPrime = circshift(LIdInwardHit,[0,0, 2^(Iter - 1)]);
                    RayArrayAlpha(LIdInwardHitPrime) = 2*pi - asin(sinAngleOfIncidence); %FIGURE OUT WHAT TO DO HERE!!!
                    RayArrayx0(LIdInwardHitPrime) = Radii(LayerId);
                    RayArrayAngDisp(LIdInwardHitPrime) = RayArrayAngDisp(LIdInwardHitPrime) + Thetaprimeprime;
                    RayArrayAmp(LIdInwardHitPrime) = ReflectionCoefficient(sinAngleOfIncidence);
                end
        %DELETING TOTALLY INTERNALLY REFLECTED RAYS (We assume that the
        %critical angle is small enough that totally internally reflected
        %rays have no chance to hit an inner shell).
        TIR = abs(imag(RayArrayAlpha)) > 0.000000000000001;
        RayArrayAlpha(TIR) = NaN;
        RayArrayx0(TIR) = NaN;
        RayArrayAmp(TIR) = NaN;
        RayArrayAngDisp(TIR) = NaN;
        IntensityArrayAngDisp(TIR) = NaN;
        IntensityArrayAmp(TIR) = NaN;
    end
end

Amp = IntensityArrayAmp(IntensityArrayAmp > 0);
AngDisp = IntensityArrayAngDisp(IntensityArrayAmp > 0);
PlaneAng = RayArrayPlaneAngle(IntensityArrayAmp > 0);
RayAng = RayArrayRayAngle(IntensityArrayAmp > 0);
[X,Y,Z] = sph2cart(AngDisp,0,1);
Xrote = X;
Yrote = cos(RayAng).*Y - sin(RayAng).*Z;
Zrote = sin(RayAng).*Y + cos(RayAng).*Z;
Xfinal = cos(PlaneAng).*Xrote -sin(PlaneAng).*Yrote;
Yfinal = sin(PlaneAng).*Xrote +cos(PlaneAng).*Yrote;
Zfinal = Zrote;
figure
axis equal
scatter3(Xfinal,Yfinal,Zfinal,'.')