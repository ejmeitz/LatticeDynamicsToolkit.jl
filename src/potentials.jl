# export StillingerWeberSilicon, LJ

# abstract type Potential end
# abstract type PairPotential <: Potential end
# abstract type ThreeBodyPotential <: Potential end

# struct LJ <: PairPotential
#     σ::Float64
#     ϵ::Float64
#     r_cut::Float64
# end


# function potential(pot::LJ, r)
#     k = (pot.σ/r)^6
#     return 4*pot.ϵ*k*(k-1)
# end


# ########################

# struct SW_Params
#     ϵ::Float64 # eV
#     σ::Float64 # Angstrom
#     a::Float64 # cutoff = a*sigma
#     λ::Float64
#     γ::Float64
#     cosθ₀::Float64
#     A::Float64
#     B::Float64
#     p::Float64
#     q::Float64
# end

# struct StillingerWeberSilicon <: ThreeBodyPotential
#     r_cut::Float64
#     params::SW_Params
#     gamma_sigma::Float64
#     lambda_epsilon::Float64
# end


# function StillingerWeberSilicon()

#     sws_params = SW_Params(2.1683, 2.0951, 1.80, 21.0, 1.20,
#                             -1/3, 7.049556277,  0.6022245584,  4.0,  0.0)
    
#     r_cut = sws_params.a*sws_params.σ
#     gamma_sigma = sws_params.σ*sws_params.γ
#     lambda_epsilon = sws_params.ϵ*sws_params.λ

#     return StillingerWeberSilicon{typeof(r_cut), typeof(gamma_sigma), typeof(lambda_epsilon)}(
#                 r_cut, sws_params, gamma_sigma, lambda_epsilon)
# end

# function pair_potential(pot::StillingerWeberSilicon, dist_ij)   
#     return Φ₂(pot.params.A, pot.params.B, pot.params.ϵ, pot.params.σ, pot.params.p, pot.params.q, pot.params.a, dist_ij)
# end

# function three_body_potential(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, dist_ij, dist_ik)
    
#     cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)

#     return Φ₃_si(pot.lambda_epsilon, pot. gamma_sigma, pot.params.cosθ₀, cosθᵢⱼₖ, pot.r_cut, dist_ij, dist_ik)

# end
 
# function Φ₂(A, B, ϵ, σ, p, q, a, rᵢⱼ)
#     return A*ϵ*(B*((σ/rᵢⱼ)^p) - ((σ/rᵢⱼ)^q)) * exp(σ/(rᵢⱼ-(a*σ)))
# end

# #Some simplifications are possible when using silicon params
# function Φ₃_si(lambda_epsilon, gamma_sigma, cosθ₀, cosθᵢⱼₖ, r_cut, rᵢⱼ, rᵢₖ)
#     return lambda_epsilon*((cosθᵢⱼₖ - cosθ₀)^2)*exp(gamma_sigma/(rᵢⱼ-r_cut))*exp(gamma_sigma/(rᵢₖ-r_cut))
# end 