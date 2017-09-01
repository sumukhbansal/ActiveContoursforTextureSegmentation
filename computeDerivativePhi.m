function dphi= computeDerivativePhi( cov_image,phi,H,Mi,Me,Xie,Xei,opt )
    
    Isize=opt.Isize;
    eps=1;
    deltaphi = deltaeps(phi,eps);
    
    k1=sum(sum(1-H));
    k3=sum(sum(H));
    
    cov_image_i=(1/(k1))*( (Mi(:)*ones(1,Isize(1)*Isize(2)))- cov_image);
    cov_image_e =(1/(k3))*( cov_image - (Me(:)*ones(1,Isize(1)*Isize(2))));
    
    
    g1 = (Mi^0.5);
    g1inv = inv(g1);
    g2 = (Me^0.5);
    g2inv = inv(g2);
    
    A1 = (g1inv' * g1inv * (-Xie) * inv(Mi))';
    vecA1 = A1(:)';
    temp_ip1 = vecA1 * cov_image_i;
    temp_dphi1 = reshape(temp_ip1,[Isize(1) Isize(2)]);
    dphi1 = deltaphi.* temp_dphi1;
    
    clear temp_ip1;
    clear temp_dphi1 ;
    clear cov_image_i;
    
    
    A2 = (g2inv' * g2inv * (-Xei) * inv(Me))';
    vecA2 = A2(:)';
    temp_ip2 = vecA2 * cov_image_e;
    temp_dphi2 = reshape(temp_ip2,[Isize(1) Isize(2)]);
    dphi2 = deltaphi .* temp_dphi2;
    
    clear temp_ip2;
    clear temp_dphi2 ;
    clear cov_image_e;
    
    dphi = dphi1 + dphi2 ;

end

