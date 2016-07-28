function output = spherical2cartesian(theta,phi,r)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x = r.*sin(theta).*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);
output = [x,y,z];
end

