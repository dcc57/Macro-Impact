%This function tells us the distance from the origin of the point source
function output=X0(L,N,n,D)
    output = (((L/N).*(n-(N/2))).^2+D^2).^(0.5);
end