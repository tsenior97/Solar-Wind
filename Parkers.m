C=-3;
f= @(r,v) v^2-log(v^2)-4*log(r)-4/r-C;
fcontour(f,[0.1 2],'MeshDensity',200)
grid on
title({'Speed of the Solar Wind'})
xlabel('r/r_c')
ylabel('v/v_c')
