using Printf
using DifferentialEquations
using LSODA
using Statistics
using Random
import Distributions: Uniform, LogNormal, fit
using PyCall
plt = pyimport("matplotlib.pyplot")
vis = PyCall.pyimport("visual") # Import python module for visualization. Module file needs to be in python lib folder.

μ_for_mean(m, σ) = log(m) - σ^2/2 # log-normal distribution

function prob_func_single(prob,i,repeat) # problem function for ODE
    prob.p["x"] = xss[i]
    prob
end

function set_pars(prob, pars) # set parameter values of a problem
    for (k, v) in pars
        prob.p[k] = v
    end
end

function prob_func(prob,i,repeat) # problem function for ensemble
    if rem(i, ncols) == 0
        prob.p["x"] = xss[ncols]
    else
        prob.p["x"] = xss[rem(i, ncols)]
    end
    prob.p["y"] = i ÷ ncols + 1
    prob
end

function transi_width(data) # computer transition width
    x = data[1,:,:,end]
    y = data[2,:,:,end]
    p5_pos = x .> minimum([0.05*maximum(x), 1.0])
    p8_pos = y .> minimum([0.05*maximum(y), 1.0])
    if maximum(x) > 300
        p5_pos .= 1
    end
    if maximum(y) > 300
        p8_pos .= 1
    end
    if all(x->x==1, p5_pos) || all(x->x==0, p5_pos) || all(x->x==1, p8_pos) || all(x->x==0, p8_pos)
        tw5, tw8 = size(data[2]), size(data[2])
    end
    p5_pos_sum = sum(p5_pos, dims=2);
    p8_pos_sum = sum(p8_pos, dims=2);
    tw5 = sum( (p5_pos_sum .< nrows) .* (p5_pos_sum .> 0) );
    tw8 = sum( (p8_pos_sum .< nrows) .* (p8_pos_sum .> 0) );
    twdp = sum( sum(p5_pos .* p8_pos, dims=2) .> 0 )
    twdn = sum( sum((p5_pos .+ p8_pos) .== 0 , dims=2) .> 0 )
    if tw5 == 0 && tw8 == 0
        xm, ym = mean(x, dims=2), mean(y, dims=2)
        if maximum(abs.(xm[2:end] - xm[1:end-1])) < (maximum(xm) - minimum(xm))*0.35
            if maximum(abs.(ym[2:end] - ym[1:end-1])) < (maximum(ym) - minimum(ym))*0.35
            tw5, tw8 = 40, 40
            end
        end
    end
    tw = max(tw5, tw8, twdp, twdn)
    return tw
end

function get_fmt_pars(pd, name_only=false) # print formated parameter values
    if name_only == true
        return rstrip(join([rpad(k,8, " ") for (k,v) in sort(collect(pd)) ], "\t"))
    else
        return rstrip(join([rpad(@sprintf("%.2f", v), 8, " ") for (k,v) in sort(collect(pd)) ], "\t"))
    end
end

mxstep(x) = maximum(abs.(mean(x, dims=2)[2:end] ./ mean(x, dims=2)[1:end-1]))
function get_mpdif(data) # compute the segregation index
    p5, p8 = data[1,:,:,end], data[2,:,:,end]
    m5, m8 = data[3,:,:,end], data[4,:,:,end]
    mpdif = (mxstep(m5) + mxstep(m8)) / (mxstep(p8) + mxstep(p5))
    return mpdif
end

include("./models/model_mmiS.jl"); # Select a model: TCR, TUR, TmiUR, TmiFB, mmiS
p_cp = copy(p); # in case we need to rerun the simulation later
prob = ODEProblem(model, u0, tspan, p) # define ODE problem
xs = LinRange(0, 100, 100) # x coordinates
ncols = 40 # number of columns in the space
xss = LinRange(xs[11], xs[81], ncols) # x coorindates for all columns
nrows = 10 # number of rows
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func) # define ensemble problem
set_pars(prob, p_cp) # in case the parameters were changed in simulations
@time sim = solve(ensemble_prob, lsoda(), reltol=1e-6, abstol=1e-6, saveat=1.0
            ,trajectories=size(xss)[1]*nrows, dt=0.01); # run simulation
num_tps = 100 # number of time points to export
d = reshape(hcat(cat(EnsembleAnalysis.componentwise_vectors_timepoint(sim,LinRange(0,1000, num_tps))..., dims=(2))...), size(u0)[1], ncols, nrows, num_tps);
println("Segregation index: ", get_mpdif(d)) # Segregation index
println("Transition width: ", transi_width(d)) # Transition width
vis.anim(d, svlist) # Visualize simulation result

