function output = theta(s,a,b,AlphaRay,x0,s1,s2,~)
%This is the angular component of the ray as a function of s
    output = acot((4 *((-1)^(s2))*a*b*x0*(a^2 - (b^2)*(x0^2))*sin(AlphaRay))/(4*(a^2)*(b^2)*(1 + exp(4*a*b*s))*(x0^2) + 4*((-1)^(s1))*a*b*exp(4*a*b*s)*x0*(a^2 + (b^2)*(x0^2))*cos(AlphaRay) + (exp(4*a*b*s) - 1)*(a^4 + (b^4)*(x0^4) + 2*(a^2)*(b^2)*(x0^2)*cos(2*AlphaRay))))- atan((1/(a^2 - (b^2)*(x0^2)))*((-1)^(s2))*(2*a*b*x0 - ((-1)^(s1))*(a^2 + (b^2)*(x0^2))*cos(AlphaRay))*csc(AlphaRay));
end