# Code taken from Wiqvist, Kalman filter for the Ornstein-Uhlenbeck model
function kalman_filter(y::Vector{Float64}, log_σ, log_c::Vector{Float64}, dt::Float64;
                       p_start=0.0)::Float64

    T = length(y)

    if log_σ isa Vector{Float64}
        log_σ = log_σ[1]
    end

    #println(loglik_est[m])
    θ_1::Float64 = exp(log_c[1])
    θ_2::Float64 = exp(log_c[2])
    θ_3::Float64 = exp(log_c[3])

    # set model
    B::Float64 = θ_1*θ_2
    A::Float64 = -θ_1
    σ::Float64 = θ_3
    C = 1
    S::Float64 = exp(log_σ)^2

    # start values for the kalman filter
    P_start = var(y) # 0.3^2
    x_hat_start = 0.0

    P_k = P_start
    x_k = x_hat_start

    loglik_est::Float64 = 0.0
    for k = 1:T
        x_k = exp(A*dt)*x_k + (1/A)*(exp(dt*A)-1)*B
        P_k = exp(A*dt)*P_k*exp(A*dt) + σ^2*(exp(2*A*dt)-1) / (2*A)

        R_k = C*P_k*C + S
        K = P_k*C*inv(R_k)
        ϵ_k = y[k]-C*x_k
        x_k = x_k + K*ϵ_k
        P_k = P_k - K*R_k*K

        loglik_est = loglik_est - 0.5*(log(det(R_k)) + ϵ_k*inv(R_k)*ϵ_k)
    end

    return loglik_est
end
