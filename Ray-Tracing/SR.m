function output = SR(R,a,b,x0,AlphaRay,s1,~,s3)
%This function tells us what value of the parameter corresponds to the ray
%intersecting a surface of radius R. May be incorrectly copied!
    output = (1/2).*a.^(-1).*b.^(-1).*log((a+(-1).*b.*R).^(-1).*(a+b.*R).^( ...
 -1).*(a.^2+b.^2.*x0.^2+2.*(-1).^s1.*a.*b.*x0.*cos(AlphaRay)).^(-1) ...
  .*(a.^4+(-1).*b.^4.*R.^2.*x0.^2+a.^2.*b.^2.*(R+(-1).*x0).*(R+ ...
  x0)+2.*(-1).^s3.*a.*b.*((R+(-1).*x0).*(R+x0).*(a.^4+(-1).*b.^4.* ...
  R.^2.*x0.^2)+(a.^2+(-1).*b.^2.*R.^2).^2.*x0.^2.*cos(AlphaRay) ...
  .^2).^(1/2)));
end