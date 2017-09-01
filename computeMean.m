function [ Me, Mi ] = computeMean( cov_image,H,opt)
% compute interior and exterior mean. 
    dim_covimage=size(cov_image);
    Isize=opt.Isize;
    
    Mi = zeros(sqrt(dim_covimage(1)));
    Me = zeros(sqrt(dim_covimage(1)));
    for i = 1:Isize(1)
        for j = 1:Isize(2)
            temp = cov_image(:,(j-1)*Isize(1) + i);
            temp1 = reshape(temp,[sqrt(dim_covimage(1)) sqrt(dim_covimage(1))]);
            Mi = Mi + (1-H(i,j))*temp1;
            Me = Me + H(i,j)*temp1;
        end
    end
    ext_no = sum(sum(H));
    Mi = Mi./(Isize(1)*Isize(2) - ext_no);
    Me = Me./(ext_no);

end

