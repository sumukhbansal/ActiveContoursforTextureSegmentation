function dlt = deltaeps(phi,eps)
  dlt = 1./(pi*eps*(1+(phi.^2/eps)));
