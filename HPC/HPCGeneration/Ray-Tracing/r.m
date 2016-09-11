function output = r(s,a,b,AlphaRay,x0,s1,~,~)
%This is the radial coordinate of the ray as a function of the parameter s
output = ((1+exp(1).^(2.*s)).^2+((-1)+exp(1).^(2.*s)).^2.*x0.^2+(-2).*((-1) ...
  +exp(1).^(4.*s)).*x0.*abs(cos(\[Alpha]))).^(-1/2).*(((-1)+exp(1).^(\
2.*s)) ...
.^2+(1+exp(1).^(2.*s)).^2.*x0.^2+(-2).*((-1)+exp(1).^(4.*s)).*x0.* ...
  abs(cos(\[Alpha]))).^(1/2);
end