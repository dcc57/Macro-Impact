clc
clear
format long
PointMatrix=[];
PointMatrix=load('20_PointsOnASphere.txt');
counter=1;
macrocounter=1;
LengthofPointmatrix=length(PointMatrix);
RayDirectionVectors=[];
SecondPhase=[];
AlphaList=[];
BetaList=[];
Flip=[]; %1 means AlphaRay->2 Pi-AlphaRay 0 means no flip
TheHolyHandArray=[]; %Angles of Incidence
Bees=[]; %Do again?
ThePlebeianThinksTheyreDone=[]; %Angle (as measured from the origin) at which the ray hits the surface of the Earth
for i=1:LengthofPointmatrix %This loop unpacks the direction vectors from a text file 'M_PointsOnASphere.txt'
    RayDirectionVectors(macrocounter,counter)=PointMatrix(i);
    if counter==3
        counter=1;
        macrocounter=macrocounter+1;
    else
        counter=counter+1;
    end
end
RE=input('What is the radius of the earth? (meters) ');
D=input('What is the distance of closest approach of the macro to the center? (meters) ');
N=input('How many points along this line should we consider? ');
Rmantle=input('What is the radius of the mantle? (meters) ');
ac=sqrt(input('Inside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is a^2? '));
bc=sqrt(input('Inside the core, assuming that the velcoity field has the form a^2-b^2*r^2, what is b^2? '));
am=sqrt(input('Outside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is a^2? '));
bm=sqrt(input('Outside the core, assuming that the velocity field has the form a^2-b^2*r^2, what is b^2? '));
vm=input('What is the velocity of p-waves just inside the mantle? ');
vc=input('What is the velocity of p-waves just outside the mantle? ');
M=LengthofPointmatrix/3;
%VARIABLE DECLARATION SECTION
%IT'S THE FINAAAAAAAL COUNTDOWN
L=2*((((RE)^2)-((D)^2))^0.5); %The length of the path through the earth.
for n=1:N-1
    x0 = X0(L,N,n,D); %X0 is the function which tells you how far from the origin the point source is, while x0 is the value for a given n
    l = ((L/N)*(n-(N/2))); %l is the distance along the trajectory from the y-axis to the point source before rotation
    Alpha = atan2(D,l); %Alpha is the angle we rotate by so that the point souce is on the x-axis
    AlphaList = cat(1,AlphaList,Alpha); %AlphaList keeps track of the Alpha's in a column vector whose length is N
    RotAlpha = [cos(Alpha),sin(Alpha),0;-sin(Alpha),cos(Alpha),0;0,0,1]; %This is the rotation matrix corresponding to a rotation about the z-axis by Alpha radians CLOCKWISE
    disp('AHHHHHHHHHHH')
        for m = 1 : M
        AlphaRayDirectionVectors=RotAlpha*[(RayDirectionVectors(m,1)+l);(RayDirectionVectors(m,2)+D);RayDirectionVectors(m,3)]; %We now rotate the ray direction vector by Alpha so that its tail lies on the x axis
        Beta=atan2(AlphaRayDirectionVectors(3),AlphaRayDirectionVectors(2)); %Beta is the angle through which we must rotate so that the direction vector lies in the xy-plane
        BetaList=cat(1,BetaList,Beta); %BetaList keeps track of all the Beta's in a column vector whose length is M
        RotBeta=[1,0,0;0,cos(Beta),sin(Beta);0,-sin(Beta),cos(Beta)]; %This is the rotation matrix corresponding to a rotation by Beta about the x-axis CLOCKWISE
        BetaAlphaRayDirectionVectors=RotBeta*[AlphaRayDirectionVectors(1);AlphaRayDirectionVectors(2);AlphaRayDirectionVectors(3)]; %These are the coordinates of the initial ray direction vectors after rotation into the xy-plane
        AlphaRay=atan2(BetaAlphaRayDirectionVectors(2),(BetaAlphaRayDirectionVectors(1)-x0)); %This is the angle that the ray direction vectors make with the x-axis
        disp('AHHHHHHHHHHH')
        if AlphaRay > pi
            AlphaRay = 2*pi - AlphaRay;
            Flip = cat(1,Flip,1);
            disp('dips')
        else
            Flip = cat(1,Flip,0);
            disp('dops')
        end
        if and(x0>Rmantle , AlphaRay > pi/2) %The ray has the potential to hit the mantle
            Smantle=SR(Rmantle,am,bm,x0,AlphaRay,1,1,0);  
            if abs(Smantle) > real(Smantle) %The ray misses the mantle
                Searth=SR(RE,am,bm,x0,AlphaRay,0,0,0);
                ThePlebeianThinksTheyreDone(m,n) = theta(Searth,am,bm,AlphaRay,x0,0,0,0);
                Bees(m,n) = 0;
            else %The ray reflects off the mantle
                TheHolyHandArray(m,n) = 2*pi - asin(GOD(Smantle,am,bm,AlphaRay,x0,1,1,0));
                Bees(m,n) = 1;
            end
        elseif and(x0>Rmantle , AlphaRay <= pi/2) %The ray does not have the potential to hit the mantle
            %(1,0,0)
            Searth=SR(RE,am,bm,x0,AlphaRay,0,0,0);
            ThePlebeianThinksTheyreDone(m,n) = theta(Searth,am,bm,AlphaRay,x0,0,0,0);
            Bees(m,n) = 0;
        elseif and(x0<=Rmantle , AlphaRay > pi/2)
            %(0,0,0)
            Bees(m,n) = 1;
        elseif and(x0<=Rmantle , AlphaRay <= pi/2) 
            %(1,0,0)
            Bees(m,n) = 1;
        end
        end
end
%{
    Note that using theta_mantle values, those are an additional rotation
    to be compensated for.  Such angle transformations are outside the
    scope of this program.
%}