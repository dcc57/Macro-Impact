function output = SX(a,b,x0,AlphaRay,s1,~,~)
% Tells us for what parameter does the ray intersect the x-axis i.e. the initial parameter value
    output = log(1 - (8*((-1)^(s1))*a*b*x0*(a^2 + (b^2)*(x0^2))*cos(AlphaRay))/(a^4 + 4*(a^2)*(b^2)*(x0^2) + (b^4)*(x0^4) + 2*a*b*x0*(2*((-1)^(s1))*(a^2 +(b^2)*(x0^2))*cos(AlphaRay) + a*b*x0*cos(2*AlphaRay))))/(4*a*b);
end

