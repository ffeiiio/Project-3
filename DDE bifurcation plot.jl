using Revise, DDEBifurcationKit, LinearAlgebra, Plots, Accessors
using BifurcationKit
const BK = BifurcationKit

# define firing_rate 
function firing_rate(x, beta)
    return 1.0 / (1.0 + exp(-beta * x))
end

# define Wilson-Cowan model
function WCvf_neuron2VF(x, xd, p)
    (; c1, c2, c3, c4, P, Q, beta) = p
    [
        -x[1] + firing_rate(c1 * xd[1][1] + c2 * xd[2][1] + P, beta),
        -x[2] + firing_rate(c3 * xd[1][2] + c4 * xd[2][2] + Q, beta)
    ]
end

function delaysF(par)
    [par.τ1, par.τ2]
end

# Set parameters
pars = (c1 = 10.0, c2 =-10.0 , c3 = 10.0, c4 = 2.0, τ1 = 5.0, τ2 = 1.0, P = -2.0, Q = -4.0, beta = 1.0)
x0 = [0.5, 3.0]

# Set P be the object 
prob = ConstantDDEBifProblem(WCvf_neuron2VF, delaysF, x0, pars, (@optic _.P))

# Use Newton iteration
optn = NewtonPar(eigsolver = DDE_DefaultEig(maxit=500))
 
opts = ContinuationPar(
    p_max = 5.0, p_min = -5.0, 
    newton_options = optn,
    ds = 0.001,              
    detect_bifurcation = 3,
    nev = 12,                
    dsmax = 0.005,          
    n_inversion = 6,         
    max_bisection_steps = 25
)


# Continuation bifurcation analyse
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = false)

# Plot bifurcation diagram
plot(br,size=(800, 700))


hpnf = BK.get_normal_form(br, 1)


brhopf = continuation(br, 1, (@optic _.τ2),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 10.0, p_min = -10.0, ds = 0.01, n_inversion = 2);
         verbosity = 0, plot = true,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)


brhopf2 = continuation(br, 2, (@optic _.τ2),
         ContinuationPar(br.contparams, detect_bifurcation = 1, dsmax = 0.01, max_steps = 100, p_max = 10.0, p_min = -10.0, ds = -0.01);
         verbosity = 2, plot = true,
         detect_codim2_bifurcation = 2,
         bothside = true,
         start_with_eigen = true)


plot(brhopf, vars = (:τ1, :τ2), legend = true)


plot(brhopf, vars = (:P, :τ1), xlims = (-2.0,-1.25), ylims = (-0,2))



plot!(brhopf2, vars = (:P, :τ1), xlims = (-3,0.0), ylims = (2.5,15))










