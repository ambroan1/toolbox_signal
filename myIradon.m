function [img,H] = myIradon( projim, thetas, f_type, f_d)

% This MATLAB function finds the inverse radon transform
% of an sinogram

% projim -projection image of size m0 x n0
% thetas - vector length n0 containing angles
% f_type - is optional and if set, it specifies the filtration to be applied during the back-projection
% f_d -  a parameter of the cut-off frequency of the filter

n = size(projim,1);  
m = size(projim,2); 
theta_rad =  (90-thetas)*pi/180;

H = NaN;

if nargin<2 || nargin==3
   error('Invalid input arguments.');

elseif nargin == 4
        H = designFilter(f_type, n, f_d);
end

[img_p,img_q] = meshgrid(linspace(n/2+1,-1*n/2,n+1),linspace(m/2+1,-1*m/2,m+1));

pmax = max(img_p,[],'all');
Pmax = sqrt(pmax^2/2);
N = floor(Pmax-2)*2;

img = zeros(N); 
[img_x,img_y] = meshgrid(linspace(-1*N/2+1,N/2,N),linspace(N/2,-1*N/2+1,1*N));

for i=1:length(thetas)  
    
    if nargin==4
        proj_fft = [projim(:,i);zeros(1,size(projim(:,i),2))]; % zero padding
        proj_filt = ifft(H.*fft(proj_fft)); % filter
        proj = proj_filt;
    else
        proj = [projim(:,i);zeros(1,size(projim(:,i),2))]; % zero padding
    end
    
    rot = img_x.*cos(theta_rad(i)) + img_y.*sin(theta_rad(i));
    new_proj = interp1(img_p(1,:),proj,rot)';
    
    img = img + new_proj;

end

img = img*pi/(2*length(thetas)); % normalization

end