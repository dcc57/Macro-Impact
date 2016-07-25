function output = r(s,a,b,AlphaRay,x0,s1,~,~)
%This is the radial coordinate of the ray as a function of the parameter s
 output = (1/(b*sqrt((a^2)*((exp(2*a*b*s) + 1)^2) + (b^2)*((exp(2*a*b*s) - 1)^2)*(x0^2) + 2*((-1)^(s1))*a*b*(exp(4*a*b*s) - 1)*x0*abs(cos(AlphaRay)))))*a*sqrt((a^2)*((exp(2*a*b*s) - 1)^2) + (b^2)*((exp(2*a*b*s) + 1)^2)*(x0^2) + 2*((-1)^(s1))*a*b*(exp(4*a*b*s) - 1)*x0*abs(cos(AlphaRay)));
end