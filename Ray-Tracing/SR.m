function output = SR(R,a,b,x0,AlphaRay,s1,s2,s3)
    output = (1/(2*a*b))log((a^4-(b^4)*(R^2)*(x0^2)+(a^2)*(b^2)*(R^2-x0^2)+((-1)^(s3))2*a*b*sqrt((R^2-x0^2)*(a^4-(b^4)*(R^2)*(x0^2))+(a^2-(b^2)*(R^2))*(x0^2)*((cos(AlphaRay))^2)))
    /((a^2-(b^2)*(R^2))((a^2+(b^2)*(R^2)+2((-1)^(s1))*a*b*x0*cos(AlphaRay)))))
end

