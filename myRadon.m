function fwrad = myRadon(img, theta)

% This MATLAB function finds the forward projection (Radon transform) of an image at specified angles
%
% img - an image matrix
% theta - vector of angles
%
% fwrad - Radon transform data (sinogram)

theta_rad =  (90-theta)*pi/180;

[img_x,img_y] = meshgrid(linspace(-1*size(img,1)/2+1/2,size(img,1)/2-1/2,size(img,1)),linspace(-1*size(img,1)/2+1/2,size(img,1)/2-1/2,size(img,1)));

xmax = max(img_x,[],'all');
ymax = max(img_y,[],'all');
Dmax = sqrt(xmax^2 +ymax^2);

[img_u,img_v] = meshgrid(-ceil(Dmax):ceil(Dmax),-ceil(Dmax):ceil(Dmax)); % [2pm+1]

fwrad = zeros(size(img_u,1),length(theta));

for i = 1:length(theta)
    img_p = img_u*cos(theta_rad(i))-img_v*sin(theta_rad(i));
    img_q = img_u*sin(theta_rad(i))+img_v*cos(theta_rad(i));
    img_tr = interp2(img_x,img_y,img,img_p,img_q,'linear',0);
    fwrad(:,i) = sum(img_tr,2);
end

end