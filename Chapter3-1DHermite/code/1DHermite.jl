function Hermite1D_Element_Matrix(hk, p=0, q=0, r=0)
    Ak = zeros(4, 4);

    Ak[1, 1] = 156*hk^2*r - 210*hk*q - 504*p;
    Ak[1, 2] = 22*hk^3*r + 42*hk^2*q - 462*hk*p;
    Ak[1, 3] = 54*hk^2*r + 210*hk*q + 504*p;
    Ak[1, 4] = -13*hk^3*r - 42*hk^2*q - 42*hk*p;

    Ak[2, 1] = 22*hk^3*r - 42*hk^2*q - 42*hk*p;
    Ak[2, 2] = 4*hk^4*r - 56*hk^2*p;
    Ak[2, 3] = 13*hk^3*r + 42*hk^2*q + 42*hk*p;
    Ak[2, 4] = -3*hk^4*r - 7*hk^3*q + 14*hk^2*p;

    Ak[3, 1] = 54*hk^2*r - 210*hk*q + 504*p;
    Ak[3, 2] = 13*hk^3*r - 42*hk^2*q + 42*hk*p;
    Ak[3, 3] = 156*hk^2*r + 210*hk*q - 504*p;
    Ak[3, 4] = -22*hk^3*r + 42*hk^2*q + 462*hk*p;

    Ak[4, 1] = -13*hk^3*r + 42*hk^2*q - 42*hk*p;
    Ak[4, 2] = -3*hk^4*r + 7*hk^3*q + 14*hk^2*p;
    Ak[4, 3] = -22*hk^3*r - 42*hk^2*q + 42*hk*p;
    Ak[4, 4] = 4*hk^4*r - 56*hk^2*p;

    return Ak;
end


function Hermite1D_Element_Integer(hk, f0, f1, f2)
    Fk = zeros(4);
    Fk[1] = 2*f0 + 4*f1;
    Fk[2] = hk * f1;
    Fk[3] = 4*f1 + 2*f2;
    Fk[4] = -Fk[2];
    Fk = Fk.*(hk^2*35.);
    return Fk;
end


function Hermite1D_ODE2_Simpson(mesh1D, f, u0, u01=NaN, p=0, q=0, r=0)
    n = length(mesh1D) - 1;
    A = zeros(2*(n+1), 2*(n+1));
    F = zeros(2*(n+1));
    u = zeros(n+1);
    u_ = zeros(n+1);

    for k = 1:n
        xk = mesh1D[k];
        hk = mesh1D[k+1] - mesh1D[k];
        Ak = Hermite1D_Element_Matrix(hk, p, q, r);
        f0, f1, f2 = f(xk), f(xk+hk/2.), f(xk+hk);
        Fk = Hermite1D_Element_Integer(hk, f0, f1, f2);
        A[2*k-1:2*k+2, 2*k-1:2*k+2] += Ak;
        F[2*k-1:2*k+2] += Fk;
    end

    if isnan(u01) == true
        A_ = A[2:end, 2:end];
        F_ = F[2:end];
        F_[1] = F_[1] - u0*A[2, 1];
        F_[2] = F_[2] - u0*A[3, 1];
        F_[3] = F_[3] - u0*A[4, 1];
        res = A\b;
        u[1] = u0;
        u_[1] = res[1];
        for i = 1:n
            u[i+1] = res[2*i];
            u_[i+1] = res[2*i+1];
        end
    else
        A_ = A[3:end, 3:end];
        F_ = F[3:end];
        F_[1] = F_[1] - u0*A[3, 1] - u01*A[3, 2];
        F_[2] = F_[2] - u0*A[4, 1] - u01*A[4, 2];
        res = A_\F_;
        u[1] = u0;
        u_[1] = u01;
        for i = 1:n
            u[i+1] = res[2*i-1];
            u_[i+1] = res[2*i];
        end
    end

    return u, u_;
end

function Hermite1D_ODE2_Interpolation(mesh1D, u, u_, p=4)
    n = length(mesh1D) - 1;
    mesh_minor = zeros(p*n + 1);
    u_minor = zeros(p*n + 1);
    u__minor = zeros(p*n + 1);
    mesh_minor[1:p:end] = mesh1D;
    u_minor[1:p:end] = u;
    u__minor[1:p:end] = u_;
    for k = 1:n
        xk = mesh1D[k];
        hk = mesh1D[k+1] - xk;
        uk00 = u[k];
        uk01 = u_[k];
        uk10 = u[k+1];
        uk11 = u_[k+1];
        for j = 1:p-1
            xi = j/p;
            mesh_minor[p*k-p+1+j] = xk + j/p*hk;
            u_minor[p*k-p+1+j] = (1-3*xi^2+2*xi^3)*uk00 + hk*(xi - 2*xi^2+xi^3)*uk01 + (3*xi^2-2*xi^3)*uk10 + hk*(-xi^2 + xi^3)*uk11;
            u__minor[p*k-p+1+j] = (-6*xi+6*xi^2)*uk00/hk + (1-4*xi+3*xi^2)*uk01 + (6*xi-6*xi^2)*uk10/hk + (-2*xi+3*xi^2)*uk11;
        end
    end
    return mesh_minor, u_minor, u__minor;
end