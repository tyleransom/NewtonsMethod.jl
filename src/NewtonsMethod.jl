module NewtonsMethod

using ForwardDiff, LinearAlgebra

function newton_root(f, f_prime; x0 = 0, tolerance=1E-7, maxiter=1_000)
    # setup the algorithm
    x_old = x0
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - f(x_old)/f_prime(x_old)
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter += 1
    end
    return (value=x_old, normdiff=normdiff, iter=iter) # A named tuple
end

function newton_root_AD(f; x0 = 0, tolerance=1E-7, maxiter=1_000)
    # setup the algorithm
    x_old = x0
    normdiff = Inf
    iter = 1
    f_prime = x -> ForwardDiff.derivative(f, x)
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - f(x_old)/f_prime(x_old)
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter += 1
    end
    return (value=x_old, normdiff=normdiff, iter=iter) # A named tuple
end

export newton_root, newton_root_AD

end
