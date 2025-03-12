using Revise, Plots
using BifurcationKit
const BK = BifurcationKit
# define function
function firing_rate(x, beta)
    return 1.0 / (1.0 + exp(-beta * x))
end

function WCvf!(dz, z, p, t = 0)
    (; c1, c2, c3, c4, P, Q, beta) = p
    u, v = z

    dz[1] = -u + firing_rate(c1 * u + c2 * v + P, beta)
    dz[2] = -v + firing_rate(c3 * u + c4 * v + Q, beta)
    dz
end
# parameter values
par_tm = (c1=10.0, c2=-10.0, c3=10.0, c4=2.0, P=-2.5, Q=-4.0, beta=1.0)

# initial condition
z0 = [0.5,3.0]

# Bifurcation Problem
prob = BifurcationProblem(WCvf!, z0, par_tm, (@optic _.P);
        record_from_solution = (x, p; k...) -> (u = x[1], v = x[2]))

# continuation options
opts_br = ContinuationPar(p_min = -10.0, p_max =1.0)

# continuation of equilibria
br = continuation(prob, PALC(tangent=Bordered()), opts_br; normC = norminf)

scene = plot(br, plotfold=false, markersize=4, legend=:topright)





