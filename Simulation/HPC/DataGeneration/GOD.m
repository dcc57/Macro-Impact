function output = GOD(s,a,b,AlphaRay,x0,s1,s2,~)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
output = 16.*(-1).^s2.*a.^3.*b.*exp(1).^(4.*a.*b.*s).*x0.*(a.^2+(-1).* ...
b.^2.*x0.^2).*(a.^2.*(1+exp(1).^(2.*a.*b.*s)).^2+b.^2.*((-1)+exp( ...
  1).^(2.*a.*b.*s)).^2.*x0.^2+2.*(-1).^s1.*a.*b.*((-1)+exp(1).^(4.* ...
  a.*b.*s)).*x0.*cos(AlphaRay)).^(-1/2).*(a.^2.*((-1)+exp(1).^(2.* ...
  a.*b.*s)).^2+b.^2.*(1+exp(1).^(2.*a.*b.*s)).^2.*x0.^2+2.*(-1) ...
.^s1.*a.*b.*((-1)+exp(1).^(4.*a.*b.*s)).*x0.*cos(AlphaRay)).^(1/2) ...
  .*(a.^4+4.*a.^2.*b.^2.*x0.^2+b.^4.*x0.^4+4.*(-1).^s1.*a.*b.*x0.*( ...
  a.^2+b.^2.*x0.^2).*cos(AlphaRay)+2.*a.^2.*b.^2.*x0.^2.*cos(2.* ...
 AlphaRay)).*(4.*a.^2.*b.^2.*(1+exp(1).^(4.*a.*b.*s)).*x0.^2+4.*( ...
  -1).^s1.*a.*b.*exp(1).^(4.*a.*b.*s).*x0.*(a.^2+b.^2.*x0.^2).*cos( ...
  AlphaRay)+((-1)+exp(1).^(4.*a.*b.*s)).*(a.^4+b.^4.*x0.^4+2.*a.^2.* ...
  b.^2.*x0.^2.*cos(2.*AlphaRay))).^(-2).*sin(AlphaRay).*(1+16.* ...
a.^2.*b.^2.*x0.^2.*(a.^2+(-1).*b.^2.*x0.^2).^2.*(4.*a.^2.*b.^2.*( ...
  1+exp(1).^(4.*a.*b.*s)).*x0.^2+4.*(-1).^s1.*a.*b.*exp(1).^(4.*a.* ...
  b.*s).*x0.*(a.^2+b.^2.*x0.^2).*cos(AlphaRay)+((-1)+exp(1).^(4.*a.* ...
  b.*s)).*(a.^4+b.^4.*x0.^4+2.*a.^2.*b.^2.*x0.^2.*cos(2.*AlphaRay))) ...
  .^(-2).*sin(AlphaRay).^2).^(-1).*(16.*a.^4.*exp(1).^(4.*a.*b.*s).* ...
  (a+(-1).*b.*x0).^2.*(a+b.*x0).^2.*(a.^2.*(1+exp(1).^(2.*a.*b.*s)) ...
  .^2+b.^2.*((-1)+exp(1).^(2.*a.*b.*s)).^2.*x0.^2+2.*(-1).^s1.*a.* ...
  b.*((-1)+exp(1).^(4.*a.*b.*s)).*x0.*cos(AlphaRay)).^(-3).*(a.^2.*( ...
  (-1)+exp(1).^(2.*a.*b.*s)).^2+b.^2.*(1+exp(1).^(2.*a.*b.*s)).^2.* ...
  x0.^2+2.*(-1).^s1.*a.*b.*((-1)+exp(1).^(4.*a.*b.*s)).*x0.*cos( ...
 AlphaRay)).^(-1).*(((-1)+exp(1).^(4.*a.*b.*s)).*(a.^2+b.^2.*x0.^2) ...
  +2.*(-1).^s1.*a.*b.*(1+exp(1).^(4.*a.*b.*s)).*x0.*cos(AlphaRay)) ...
  .^2+256.*a.^6.*b.^2.*exp(1).^(8.*a.*b.*s).*x0.^2.*(a.^2+(-1).* ...
 b.^2.*x0.^2).^2.*(a.^2.*(1+exp(1).^(2.*a.*b.*s)).^2+b.^2.*((-1)+ ...
  exp(1).^(2.*a.*b.*s)).^2.*x0.^2+2.*(-1).^s1.*a.*b.*((-1)+exp(1).^( ...
  4.*a.*b.*s)).*x0.*cos(AlphaRay)).^(-1).*(a.^2.*((-1)+exp(1).^(2.* ...
  a.*b.*s)).^2+b.^2.*(1+exp(1).^(2.*a.*b.*s)).^2.*x0.^2+2.*(-1) ...
.^s1.*a.*b.*((-1)+exp(1).^(4.*a.*b.*s)).*x0.*cos(AlphaRay)).*( ...
a.^4+4.*a.^2.*b.^2.*x0.^2+b.^4.*x0.^4+4.*(-1).^s1.*a.*b.*x0.*( ...
a.^2+b.^2.*x0.^2).*cos(AlphaRay)+2.*a.^2.*b.^2.*x0.^2.*cos(2.* ...
AlphaRay)).^2.*(4.*a.^2.*b.^2.*(1+exp(1).^(4.*a.*b.*s)).*x0.^2+4.* ...
  (-1).^s1.*a.*b.*exp(1).^(4.*a.*b.*s).*x0.*(a.^2+b.^2.*x0.^2).*cos( ...
  AlphaRay)+((-1)+exp(1).^(4.*a.*b.*s)).*(a.^4+b.^4.*x0.^4+2.*a.^2.* ...
  b.^2.*x0.^2.*cos(2.*AlphaRay))).^(-4).*sin(AlphaRay).^2.*(1+16.* ...
  a.^2.*b.^2.*x0.^2.*(a.^2+(-1).*b.^2.*x0.^2).^2.*(4.*a.^2.*b.^2.*( ...
  1+exp(1).^(4.*a.*b.*s)).*x0.^2+4.*(-1).^s1.*a.*b.*exp(1).^(4.*a.* ...
  b.*s).*x0.*(a.^2+b.^2.*x0.^2).*cos(AlphaRay)+((-1)+exp(1).^(4.*a.* ...
  b.*s)).*(a.^4+b.^4.*x0.^4+2.*a.^2.*b.^2.*x0.^2.*cos(2.*AlphaRay))) ...
  .^(-2).*sin(AlphaRay).^2).^(-2)).^(-1/2);
end

