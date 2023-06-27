% inspired by Abinit subroutine coeffs_gausslegint
function [omega_mesh, mesh01, weights] = imaginary_omega_mesh(number)
    omega_mesh = zeros(1, number);
    [mesh01, weights] = gauss_legendre_coeffs(0.0, 1.0, number); % Replace $ \int_0^\infty dx f(x) $ with $ \int_0^1 dz f(1/z - 1)/z^2 $.
    for index = 1:number
        omega_mesh(index) = 1.0/mesh01(index) - 1.0;
    end
end

function [mesh, weights] = gauss_legendre_coeffs(minVal, maxVal, number)
    tolerance = 1e-13;
    length = (maxVal - minVal) / 2;
    mean = (maxVal + minVal) / 2;
    mesh = zeros(1, number);
    weights = zeros(1, number);
    for index = 1:floor((number + 1) / 2)
        z = cos(pi*(index - 0.25)/(number + 0.5));
        flag = 1;
        while flag == 1
            p1 = 1.0;
            p2 = 0.0;
            for j = 1:number
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3) / j;
            end
            pp = number*(p2 - z*p1) / (1.0 - z*z);
            z1 = z;
            z = z1 - p1/pp;
            flag = (abs(z - z1) > tolerance);
        end
        mesh(index) = mean - length*z;
        mesh(number + 1 - index) = mean + length*z;
        weights(index) = 2.0*length / ((1.0 - z*z) * (pp*pp));
        weights(number + 1 - index) = weights(index);
    end
end