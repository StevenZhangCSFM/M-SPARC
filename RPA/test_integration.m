% test integration
number = 12;
[omega_mesh, mesh01, weights] = imaginary_omega_mesh(number);
values = exp(-omega_mesh);
integral = sum(values./(mesh01.^2) .* weights);