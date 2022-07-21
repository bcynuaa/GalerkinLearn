function Galerkin1D_ODE1_Simpson(mesh1D, u0, f)
    n = length(mesh1D) - 1;
    u = zeros(n+1);
    u[1] = u0;
    h0 = mesh1D[2] - mesh1D[1];
    u[2] = u0 + h0/3. * ( f(mesh1D[1]) + 2. * f(mesh1D[1] + h0/2.));
    for i = 2:n-1
        hi0 = mesh1D[i] - mesh1D[i-1];
        hi1 = mesh1D[i+1] - mesh1D[i];
        temp = 2*hi0/3 * f(mesh1D[i] - hi0/2.);
        temp += (hi0+hi1)/3 * f(mesh1D[i]);
        temp += 2*hi1/3 * f(mesh1D[i] + hi1/2.);
        u[i+1] = u[i-1] + temp;
    end
    hn = mesh1D[n+1] - mesh1D[n];
    u[n+1] = u[n] + hn/3. * ( f(mesh1D[n]) + 2. * f(mesh1D[n] + hn/2.) );
    return u;
end