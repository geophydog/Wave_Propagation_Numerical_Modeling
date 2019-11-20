function [Ke, Me] = KeMe( dx, dz, c )
gamma = dz / dx;
Ke = zeros(4, 4);
Me = zeros(4, 4);
for i = 1: 4
    Ke(i, i) = gamma + 1. / gamma;
    Me(i, i) = 1.;
end
Ke(1, 2) = gamma/2. - 1./gamma; Ke(1, 3) = -0.5*(gamma+1./gamma);
Ke(1, 4) = 1./(2.*gamma)-gamma;
Ke(2, 1) = Ke(1, 2); Ke(2, 3) = 1./(2*gamma)-gamma; Ke(2, 4) = -0.5*(gamma+1./gamma);
Ke(3, 1) = Ke(1, 3); Ke(3, 2) = Ke(2, 3); Ke(3, 4) = gamma/2.-1./gamma;
Ke(4, 1) = Ke(1, 4); Ke(4, 2) = Ke(2, 4); Ke(4,3) = Ke(3, 4);

Ke = Ke * (c^2/3.);
Me = Me * (dx*dz/4.);
end

