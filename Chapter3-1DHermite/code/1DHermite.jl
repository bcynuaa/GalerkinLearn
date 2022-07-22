function Hermite1D_Element_Matrix(hk, p=0, q=0, r=0)
    Ak = zeros(4, 4);

    Ak[1][1] = 156*hk^2*r - 210*hk*q - 504*p;
    Ak[1][2] = 22*hk^3*r + 42*hk^2*q - 462*hk*p;
    Ak[1][3] = 54*hk^2*r + 210*hk*q + 504*p;
    Ak[1][4] = -13*hk^3*r - 42*hk^2*q - 42*hk*p;

    Ak[2][1] = 22*hk^3*r - 42*hk^2*q - 42*hk*p;
    Ak[2][2] = 4*hk^4*r - 56*hk^2*p;
    Ak[2][3] = 13*hk^3*r + 42*hk^2*q + 42*hk*p;
    Ak[2][4] = -3*hk^4*r - 7*hk^3*q + 14*hk^2*p;

    Ak[3][1] = 54*hk^2*r - 210*hk*q + 504*p;
    Ak[3][2] = 13*hk^3*r - 42*hk^2*q + 42*hk*p;
    Ak[3][3] = 156*hk^2*r + 210*hk*q - 504*p;
    Ak[3][4] = -22*hk^3*r + 42*hk^2*q + 462*hk*p;

    Ak[4][1] = -13*hk^3*r + 42*hk^2*q - 42*hk*p;
    Ak[4][2] = -3*hk^4*r + 7*hk^3*q + 14*hk^2*p;
    Ak[4][3] = -22*hk^3*r - 42*hk^2*q + 42*hk*p;
    Ak[4][4] = 4*hk^4*r - 56*hk^2*p;

    return Ak;
end

function Hermite1D_Element_Integer(xk, hk, f)
    Fk = zeros(4);
    Fk[1] = 2*f(xk) + 4*f(xk + hk/2);
    Fk[2] = hk * f(xk + hk/2);
    Fk[3] = 4*f(xk + hk/2) + 2*f(xk + hk);
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
        Fk = Hermite1D_Element_Integer(xk, hk, f);
        A[2*k-1:2*k+2] += Ak;
        F[2*k+1-1:2*k+2] += Fk;
    end

    if isnan(u01) == True
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
        res = A\b;
        u[1] = u0;
        u_[1] = u01;
        for i = 1:n
            u[i+1] = res[2*i-1];
            u_[i+1] = res[2*i];
        end
    end

    return u, u_;
end