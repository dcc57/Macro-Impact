clc
clear
format long
%EvenlyDistributedPoints=[];
%EvenlyDistributedPoints=load('130_Points_On_A_Sphere.txt');
counter=1;
macrocounter=1;
%LengthofPointmatrix=length(EvenlyDistributedPoints);
%RayDirectionVectors=zeros(LengthofPointmatrix/3,3);
SecondPhase=[];
%for i = 1 : LengthofPointmatrix %This loop unpacks the direction vectors from a text file 'M_PointsOnASphere.txt'
%    RayDirectionVectors(macrocounter,counter) = EvenlyDistributedPoints(i);
%    if counter == 3
%        counter = 1;
%        macrocounter = macrocounter + 1;
%    else
%        counter = counter + 1;
%    end
%end
RE = 6371; %input('What is the radius of the earth? (meters) ');
D = 3000; %input('What is the distance of closest approach of the macro to the center? (meters) ');
N = 100; %input('How many points along this line should we consider? ') + 1;
M = 100; %input('How many rays will each point emit? ');
Rmantle = 3480; %input('What is the radius of the mantle? (meters) ');
ac = 11.1502; %sqrt(input('Inside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is a^2? '));
bc = 2.69692*(10^(-7)); %sqrt(input('Inside the core, assuming that the velcoity field has the form a^2-b^2*r^2, what is b^2? '));
am = 15.4938; %sqrt(input('Outside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is a^2? '));
bm = 1.40336*(10^(-7)); %sqrt(input('Outside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is b^2? '));
RayDirectionVectors = zeros(M,N,3);
for n = 1 : N
    for m = 1 : M
        Th = 2*pi*rand;
        Ph = asin(-1+2*rand);
        [x,y,z] = sph2cart(Th,Ph,1);
        RayDirectionVectors(m,n,:) = [x,y,z];
    end
end
vm = (am^2)-(bm^2)*(Rmantle)^2;
vc = (ac^2)-(bc^2)*(Rmantle)^2;
CriticalAngle = asin(vc/vm);
%M = LengthofPointmatrix/3; %The number of rays propagating from each point
LaunchAngles = zeros(M,N); %Angles of Incidence
AngleAtCoreMantleBoundary = zeros(M,N); %Angular displacement for rays incident on the core-mantle boundary
AngleAtSurface=zeros(M,N); %Angular displacement for rays incident on the surface of the earth
HitCore = zeros(M,N); %Do again?
HitSurface = zeros(M,N); %Did the ray get totally internally reflected?
AlphaList = zeros(1,N); %A running list of Alpha's. These are the angles by which the point sources are rotated.
BetaList = zeros(M,N); %A running list of Beta's.
ImpactPositionList = zeros(M*N,3); %The final xyz coordinates of the rays. This is an MxN array with a triplet for each entry.
L = 2*((((RE)^2)-((D)^2))^0.5); %The length of the path through the earth.
for n=1:N-1
    x0 = X0(L,N,n,D); %X0 is the function which tells you how far from the origin the point source is, while x0 is the value for a given n
    l = ((L/N)*(n-(N/2))); %l is the distance along the trajectory from the y-axis to the point source before rotation
    Alpha = atan2(D,l); %Alpha is the angle we rotate by so that the point souce is on the x-axis
    AlphaList(n) = Alpha; %AlphaList keeps track of the Alpha's in a row vector whose length is N
    RotAlpha = [cos(Alpha),sin(Alpha),0;-sin(Alpha),cos(Alpha),0;0,0,1]; %This is the rotation matrix corresponding to a rotation about the z-axis by Alpha radians CLOCKWISE
        for m = 1 : M
        AlphaRayDirectionVectors(m,n,:) = RotAlpha*[(RayDirectionVectors(m,n,1)+l);(RayDirectionVectors(m,n,2)+D);RayDirectionVectors(m,n,3)]; %We now rotate the ray direction vector by Alpha so that its tail lies on the x axis
        Beta = atan2(AlphaRayDirectionVectors(m,n,3),AlphaRayDirectionVectors(m,n,2)); %Beta is the angle through which we must rotate so that the direction vector lies in the xy-plane
        BetaList(m,n) = Beta; %BetaList keeps track of all the Beta's in a column vector whose length is M
        RotBeta = [1,0,0;0,cos(Beta),sin(Beta);0,-sin(Beta),cos(Beta)]; %This is the rotation matrix corresponding to a rotation by Beta about the x-axis CLOCKWISE
        BetaAlphaRayDirectionVectors(m,n,:) = RotBeta*[AlphaRayDirectionVectors(m,n,1);AlphaRayDirectionVectors(m,n,2);AlphaRayDirectionVectors(m,n,3)]; %These are the coordinates of the initial ray direction vectors after rotation into the xy-plane
        AlphaRay = atan2(BetaAlphaRayDirectionVectors(m,n,2),(BetaAlphaRayDirectionVectors(m,n,1)-x0)); %This is the angle that the ray direction vectors make with the x-axis
        if and(x0 > Rmantle , AlphaRay > pi/2) %The ray has the potential to hit the mantle
            AlphaRay = pi - AlphaRay;
            Smantle = SR(Rmantle,am,bm,x0,AlphaRay,0,0,1);
            if abs(imag(Smantle))>0 %The ray misses the mantle
                Searth=SR(RE,am,bm,x0,AlphaRay,0,0,0);
                AngleAtSurface(m,n) = theta(Searth,am,bm,AlphaRay,x0,0,0,0);
                HitSurface(m,n) = 1;
            else %The ray reflects off the mantle
                LaunchAngles(m,n) = 2*pi - asin(GOD(Smantle,am,bm,AlphaRay,x0,0,0,1));
                AngleAtCoreMantleBoundary(m,n) = theta(Smantle,am,bm,AlphaRay,x0,0,0,1);
                HitCore(m,n) = 1;
            end
        elseif and(x0>Rmantle , AlphaRay <= pi/2) %The ray does not have the potential to hit the mantle
            Searth=SR(RE,am,bm,x0,AlphaRay,1,0,0);
            AngleAtSurface(m,n) = theta(Searth,am,bm,AlphaRay,x0,1,0,0);
            HitSurface(m,n) = 1;
        elseif and(x0<=Rmantle , AlphaRay > pi/2)
            AlphaRay = pi - AlphaRay;
            Smantle = SR(Rmantle,ac,bc,x0,AlphaRay,0,0,0);
            god = GOD(Smantle,ac,bc,AlphaRay,x0,0,0,0);
            angleofincidence = asin(god);
            if CriticalAngle > angleofincidence
                LaunchAngles(m,n) = asin((vm/vc)*god);
                AngleAtCoreMantleBoundary(m,n) = theta(Smantle,ac,bc,AlphaRay,x0,0,0,0);
                HitCore(m,n) = 1;
            end
        elseif and(x0<=Rmantle , AlphaRay <= pi/2) 
            Smantle = SR(Rmantle,ac,bc,x0,AlphaRay,1,0,0);
            god = GOD(Smantle,ac,bc,AlphaRay,x0,1,0,0);
            angleofincidence = asin(god);
            if CriticalAngle > angleofincidence
                LaunchAngles(m,n) = asin((vm/vc)*god);
                AngleAtCoreMantleBoundary(m,n) = theta(Smantle,ac,bc,AlphaRay,x0,1,0,0);
                HitCore(m,n) = 1;
            end
        end
        end
end
for n = 1 : N - 1
    for m = 1 : M
        if HitCore(m,n) == 1 %note that the angle of departure is always between 0 and pi/2
            Searth = SR(RE,am,bm,Rmantle,LaunchAngles(m,n),1,0,0);
            AngleAtSurface(m,n) = theta(Searth,am,bm,LaunchAngles(m,n),Rmantle,1,0,0)+AngleAtCoreMantleBoundary(m,n);
            HitSurface(m,n) = 1;
        end
    end
end
NumberOfImpacts = 0;
for n = 1 : N - 1
    for m = 1 : M
        if HitSurface(m,n) == 1
            NumberOfImpacts = NumberOfImpacts + 1;
            ImpactPositionBeforeRotation = [cos(AngleAtSurface(m,n));sin(AngleAtSurface(m,n));0];
            RotBackBeta = [1,0,0;0,cos(BetaList(m,n)),-sin(BetaList(m,n));0,sin(BetaList(m,n)),cos(BetaList(m,n))];
            ImpactPositionAfterFirstRotation = RotBackBeta * ImpactPositionBeforeRotation;
            RotBackAlpha = [cos(AlphaList(n)),-sin(AlphaList(n)),0;sin(AlphaList(n)),cos(AlphaList(n)),0;0,0,1];
            ImpactPosition = RotBackAlpha * ImpactPositionAfterFirstRotation;
            for i = 1 : 3
                ImpactPositionList(NumberOfImpacts,i) = ImpactPosition(i);
            end
        end
    end
end
ImpactPositionList = ImpactPositionList(1:NumberOfImpacts , : );
for k = 1 : NumberOfImpacts
    [TH,PH,R]=cart2sph(ImpactPositionList(k,1),ImpactPositionList(k,2),ImpactPositionList(k,3));   
    ImpactPositionListPolar(k,:) = [TH,PH,R];
end

%Now we need to convert this into a pretty picture. This section is slow,
%and could use some work...
figure
scatter3(ImpactPositionList(:,1),ImpactPositionList(:,2),ImpactPositionList(:,3),'.')
PolarDensity = 60;
AzimuthalDensity = 60;
NumberOfSamplePoints = (PolarDensity) * (AzimuthalDensity - 1);
SampleAngle = pi/30;
PointsToBeColored = zeros(AzimuthalDensity * PolarDensity,3);
CartPointsToBeColored = zeros((AzimuthalDensity - 1) * PolarDensity,3);
RotatedPositionList = zeros(NumberOfImpacts,3);
NumberOfRaysIncidentInDisk = zeros(NumberOfSamplePoints,1);
k = 0;
for i = 1 : AzimuthalDensity - 1
    for j = 1 : PolarDensity
        k = k + 1;
        PointsToBeColored(k,:) = [((pi)/(AzimuthalDensity - 1)) * i ,((2*pi)/(PolarDensity))* j, 1];
        [x,y,z] = sph2cart(PointsToBeColored(k,1),PointsToBeColored(k,2),PointsToBeColored(k,3));
        CartPointsToBeColored(k,:) = [x,y,z];
    end
end
for i = 1 : NumberOfSamplePoints
    PHI = atan2((sqrt((CartPointsToBeColored(i,1))^2+(CartPointsToBeColored(i,2))^2)),(CartPointsToBeColored(i,3)));
    for m = 1 : NumberOfImpacts
        RotatedPositionList(m,:) = T(-CartPointsToBeColored(i,2),CartPointsToBeColored(i,1),0,-PHI)*transpose(ImpactPositionList(m,:));
    end
    for n = 1 : NumberOfImpacts
        phi = atan2((sqrt((RotatedPositionList(n,1))^2+(RotatedPositionList(n,2))^2)),(RotatedPositionList(n,3)));
        RotatedPositionListPolar(:) = [Th,Ph,R];
        if abs(phi) <= SampleAngle
            NumberOfRaysIncidentInDisk(i,1) = NumberOfRaysIncidentInDisk(i,1) + 1;
        end
    end
end
figure
scatter3(CartPointsToBeColored(:,1),CartPointsToBeColored(:,2),CartPointsToBeColored(:,3),200,log(NumberOfRaysIncidentInDisk(:)),'filled')