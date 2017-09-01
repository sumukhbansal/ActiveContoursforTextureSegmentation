function cov_image = computeMVI( I,opt )
% Compute manifold valued image by computing covarience of the patch around each pixel. 
Isize=opt.Isize;
dim_patch=opt.dim_patch;
s=floor(dim_patch/2);
cov_image = zeros((2*s+1)^4,Isize(1)*Isize(2));
I_padded=padarray(double(I),[s s],'symmetric','both');
% I_pad=padarray(double(I),[1 1],'symmetric','both');
for i = 1:Isize(1)
    for j = 1:Isize(2)
        temp = I_padded(i:i+2*s,j:j+2*s);
        temp1 = temp(:);
        tcov = temp1 * temp1';
        tcov1 = tcov(:);
        cov_image(:,(j-1)*Isize(1) + i) = tcov1;
    end
end
end

