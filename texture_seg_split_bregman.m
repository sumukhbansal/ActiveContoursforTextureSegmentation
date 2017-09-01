function [I,phi]= texture_seg_split_bregman( I,opt)
% Main function for texture segmenatation using split bregman technique.
% Inputs: 
% I             : Input image
% opt.Isize     : image size
% opt.lambda    : parameter for split-bregman algorithm
% opt.mu        : parameter for split-bregman algorithm  
% opt.dim_patch : patch dimentions
% opt.max_itr   : maximum iterations
% Output        : Resized image and Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resize and convert color image to gray image
Isize =opt.Isize;
if (length(size(I))==3)
    I = imresize(rgb2gray(I),Isize);
else
    I = imresize(I,Isize);
end

lamda=opt.lambda;
mu=opt.mu;

% default parameters and initialize phi
u = ones(Isize); 
u(30:Isize(1)-30,30:Isize(2)-30) = 0;
dx = zeros(Isize);
bx = zeros(Isize);
dy = zeros(Isize);
by = zeros(Isize);
eps =1;

% compute covarience of the patch around each pixel 
cov_image = computeMVI( I,opt );

phi = u;
err=1;
dist = 0;
itr=0;

% main loop
while(err>0 && itr< opt.max_itr)
    
    H = Hstep(phi,eps);
    [ Me, Mi ] = computeMean( cov_image,H,opt);
    Xie = log_pd(Mi,Me);
    Xei = log_pd(Me,Mi);
    
    dist_old = dist;
    dist = sqrt(ippdn(Xie,Xie,Mi));
    err=dist-dist_old;
    disp(err);
    if (err<0)
        return
    end
    dphi= computeDerivativePhi( cov_image,phi,H,Mi,Me,Xie,Xei,opt );
    
    alpha = circshift(dx,[0,1]) - dx - circshift(bx,[0,1]) + bx + circshift(dy,[1 ,0]) - dy - circshift(by,[1,0]) + by;
    beta  =(1/4)*( circshift(u,[0,1]) + circshift(u,[0,-1]) + circshift(u,[1,0]) + circshift(u,[-1,0]) + alpha - (lamda/mu)* dphi);
    
    [ux, uy]=gradient(u);
    dx = sign(ux + bx)* max( ((ux + bx).^2).^(.5) - (1/lamda), zeros(Isize));
    dy = sign(uy + by)* max( ((uy + by).^2).^(.5) - (1/lamda), zeros(Isize));
    bx = bx + ux - dx;
    by = by + uy - dy;
    u = max( min(beta,ones(Isize)),zeros(Isize));
    phi=u;
    itr = itr+1;
    imshow(I);hold on;
    [C,h] = contour(phi,'r','LineWidth',1);
    pause(.1);
    
end

end

