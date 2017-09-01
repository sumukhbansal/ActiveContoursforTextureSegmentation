function ip = ippdn(Xp,Yp,p)
g = (p^0.5);
ginv = inv(g);
ip = trace(ginv * Xp * inv(p) * Yp * ginv');
