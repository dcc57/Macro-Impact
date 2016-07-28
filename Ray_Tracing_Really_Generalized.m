


clc
clear
format long

%DEFINITION OF PARAMETERS

D = 3000; %input('What is the distance of closest approach of the macro to the center? (meters) ');
N = 100; %input('How many points along this line should we consider? ') + 1;
M = 100; %input('How many rays will each point emit? ');
Layers = 3; %input('How many discrete strata is the Earth composed of? ');
Radii = [0,1221,3480,6371]; %What are the radii of the strata? (Layers + 1)
a = [11.2633,10.6835,15.4938].^(1/2); %first entry corresponds to the core
b = [1.55946*10^(-7),2.11853*10^(-7),1.40336*10^(-7)].^(1/2); %first entry corresponds to the core
WeaknessOfNeglectedRays = 10^(-1); %input('What fraction of the initial ray energy must a ray contain to not be thrown out? ');
AbsoluteCutoff = 5; %Number of times a ray can bounce before it's thrown out.
RayDirectionVectors = zeros(M,N,3);
AlphaRay = zeros(M,N);
for n = 1 : N
    for m = 1 : M
        Th = 2*pi*rand;
        Ph = asin(-1+2*rand);
        [x,y,z] = sph2cart(Th,Ph,1);
        RayDirectionVectors(m,n,:) = [x,y,z];
    end
end

vBoundaries = zeros(layers-1,2); %v(:,1) corresponds to the inner velocity on the :th layer, and v(:,2) corresponds to the outer velocity on the :th layer
ReflectionCoefficients = zeros(layers - 1,2); %these correspond to 4 functions of the angle of incidence
TransmissionCoefficients = zeros(layers - 1,2); %these correspond to 4 functions of the angle of incidence
for i = 1 : Layers - 1
    vBoundaries(i,1) = (a(i)^2)-(b(i)^2)*(Radii(i))^2; %inner
    vBoundaries(i,2) = (a(i+1)^2)-(b(i+1)^2)*(Radii(i))^2; %outer
end

%RAY INITIAL TRAJECTORIES

AlphaList = zeros(1,N); %A running list of Alpha's. These are the angles by which the point sources are rotated.
BetaList = zeros(M,N); %A running list of Beta's.
L = 2*((((RE)^2)-((D)^2))^0.5); %The length of the path through the earth.
%Each initial ray can correspond to (assuming each ray hits a boundary)
%2^(AbsoluteCutoff-1) rays after AbsoluteCutoff iterations.
x0 = zeros(N,1);
x0(:,1) = X0(L,(N+1),:,D); %X0 is the function which tells you how far from the origin the point source is, while x0 is the value for a given n
counter = 0;
for n=1:N
    l = ((L/(N+1))*(n-((N+1)/2))); %l is the distance along the trajectory from the y-axis to the point source before rotation
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
    end
end
%We now have the directions of the rays for the first round of propagation.
%Note that the children of each ray will inherit the parent ray's AlphaList
%and BetaList, so there is no additional complication from that. There are
%two ways a ray can propagate within a given strata. First, it can
%propagate to the next layer. Second, it can propagate to the previous
%layer. Thus, for each ray, we have only two cases to consider.

%RAY PROPAGATION

NumberOfRays = zeros(AbsoluteCutoff,1);
NumberOfRays(1,1) = N*M;

RayArray = zeros(M,N,2^(AbsoluteCutoff - 1),3);
RayArray(1 : M,1 : N,1,1) = X0(L,(N+1),:,D);
RayArray(1 : M,1 : N,1,2) = AlphaRay(:,:);
RayArray(1 : M,1 : N,1,3) = 1;
for IterationCounter = 1 : AbsoluteCutoff
    [1Id, 2Id, 3Id, 4Id] = find( RayArray );
    s2in(:,:,:) = RayArray(:,:,3Id,4Id);
    s2in(imag(s2in) == 0) = 0;
    s2out(:,:,:) = 
end


for IterationCounter = 1 : 5
    for LayerCounter = 2 : Layers + 1
        for n = 1 : NumberOfRays(IterationCounter,1)
            if and(x0(n) < Radii(LayerCounter),x0(n) >= Radii(LayerCounter - 1))
                if %angle pointing outwards - only one case
                    if LayerCounter == Layers + 1
                        %Store Ray in Final Data Spot, and leave out of
                        %computations!
                        NumberOfRays(IterationCounter + 1,1) = NumberOfRays(IterationCounter,1) - 1;
                    elseif IterationCounter < AbsoluteCutoff
                        NumberOfRays(IterationCounter + 1,1) = NumberOfRays(IterationCounter,1) + 1;
                    end
                elseif %angle pointing inwards - two cases
                    if %s(lower layer) is real
                    else %s(lower layer) is complex
                        if LayerCounter == Layers + 1
                            %Store Ray in Final Data Spot, and leave out of
                            %computations!
                            NumberOfRays(IterationCounter + 1,1) = NumberOfRays(IterationCounter,1) - 1;
                        elseif IterationCounter < AbsoluteCutoff
                            NumberOfRays(IterationCounter + 1,1) = NumberOfRays(IterationCounter,1) + 1;
                        end
                    end
                end
            end
        end
    end
end












