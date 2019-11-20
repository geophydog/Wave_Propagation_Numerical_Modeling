% No. of time steps.
nt = 201;
% Temporal interval.
dt = 0.001;
% Spatial interval on X axis.
dx = 1.;
% Spatial interval on Z axis.
dz = 1.;
% No. of spatial steps on X axis.
nx = 111;
% No. of spatial steps on Z axis.
nz = 111;
N = floor(nz * nz);
% Position of source
sx = floor(nx/2);
sz = floor(nz/2);
x = (0: nx-1) * dx;
z = (0: nz-1) * dz;
% Acoustic velocity in m/s.
c = 500.;
% Center frequency of source Ricker wavelet.
f0 = 30;
t = (0: nt) * dt;
% Source time function: Ricker wavelet.
t0 = 35*dt;
s = (1-2*(pi*f0*(t-t0)).^2) .* exp(-(pi*f0*(t-t0)).^2);

% Position of source.
f = zeros(nz, nz);
f(sz, sx) = 1.;
u = zeros(nt, N);
f = reshape(f', [N, 1]);

% Elementary matrices and global matrices.
[Ke,Me] = KeMe(dx, dz, c);
[K, M] = globalKM(nz, Ke, Me);

% Invert global mass matrix.
%INVM = inv(M);
for i = 1: nz
    M(i, i) = 1. / M(i, i);
end

% Iterations for finding solutions.
for i = 2: (nt-1)
    u(i+1, :) = 2.*u(i, :)' - u(i-1, :)' + dt^2 * M *...
        (-K * u(i, :)' + f*s(i));
end

% Plot.
ti = 20;
tt = floor(nt/ti);
figure(1);
for i = 1: ti
    imagesc(x, z, reshape(u(floor(i*tt), :), [nz, nz])');
    hold on;
    scatter(sx*dx, sz*dz, 'filled', 'o');
    hold on;
    xlabel('X [m]');
    ylabel('Z [m]');
    set(gca, 'fontsize', 20, 'fontweight', 'bold', 'ydir', 'normal');
end