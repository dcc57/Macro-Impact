clear

%GENERATING THE RANDOM POINTS ON A SPHERE

M = 7;
N = 5;
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

D = 2000;
Radii = [0,1221,3480,6371];
AbsoluteCutoff = 5;

Layers = length(Radii);

L = 2*((((Radii(Layers))^2)-((D)^2))^0.5);
l = (I-((M+1)/2)).*(L/(M+1));
RayArray = zeros(M,N,2^(AbsoluteCutoff - 1),3); % = (point M, ray N, daughter P, (x,AlphaRay,Amplitude))
Alpha = atan2(D,transpose(l));
for i = 1 : M
    RotMatrixAlpha = [cos(Alpha(i)),sin(Alpha(i)),0;-sin(Alpha(i)),cos(Alpha(i)),0;0,0,1];
    for j = 1 : N
        RotRand = RotMatrixAlpha*[RandX(i,j)+l(i);RandY(i,j)+D;RandZ(i,j)];
        Beta(i,j) = atan2(RotRand(3),RotRand(2));
        RotMatrixBeta = [1,0,0;0,cos(Beta(i,j)),sin(Beta(i,j));0,-sin(Beta(i,j)),cos(Beta(i,j))];
        Vec = RotMatrixBeta*RotRand;
        AlphaRay(i,j) = atan2(Vec(2),Vec(1));
    end
end
x0 = transpose(repmat(X0(L,(M+1),I,D), [N,1]));
RayArray(:,:,1,1) = x0(:,:);
RayArray(:,:,1,2) = AlphaRay(:,:);
RayArray(:,:,1,3) = 1;

%PROPAGATION

LeftRays = zeros(M,N,2^(AbsoluteCutoff - 1),1);
LayerRays = zeros(M,N,2^(AbsoluteCutoff - 1),Layers - 1);

for Iteration = 1 : AbsoluteCutoff
    LeftRays = ( RayArray(:,:,:,2) > pi / 2 );
    RayArray(:,:,:,LeftRays) = pi - RayArray(:,:,:,LeftRays);
    for LayerCounter = 1 : Layers - 1;
        LayerRays(:,:,:,LayerCounter) = ( and( RayArray(:,:,:,1) > Radii(LayerCounter), RayArray(:,:,:,1) < Radii(LayerCounter + 1) ) );
    end
    RayArray(:,:,:,2)
end









