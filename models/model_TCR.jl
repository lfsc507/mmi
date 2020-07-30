function model(du,u,p,t)
    ts = [50.0, 100.0, 200.0, 900.0]
    freq = 1.0 * 0.2
    amp = p["amp"]
    if t < ts[1]
        p["maxA"] = 0
    elseif t < ts[3]
        p["maxA"] = 2.3
    else
        p["maxA"] = 2.3*0.65
        p["gA"] = 0.005
    end
    if t < ts[2]
        p["maxF"] = 0
    elseif t < ts[3]
        p["maxF"] = 2.1
    else
        p["maxF"] = 2.1*0.65
        p["gF"] = 0.005
    end
    if t < ts[4]
        p["muA"] = 0
        p["muF"] = 0
    else
        Random.seed!(floor(Int, (t-ts[4])*freq)*p["y"]*floor(Int, p["x"]))
        p["muA"] =   sign(sin((t-ts[4])*pi*freq)) * amp * randn()*u[5]
        p["muF"] =  sign(sin((t-ts[4])*pi*freq)) * amp * randn()*u[6]
    end
    P5, P8, R5, R8, A, F, r2, c15, c25, c35, r9, c18, c28, c38 = u
    x, gA, gF, maxA, maxF, s05, s08, K5A, K8F, n5A, n8F, K58, n58, K85, n85  = p["x"], p["gA"], p["gF"], p["maxA"], p["maxF"], p["s05"], p["s08"], p["K5A"], p["K8F"], p["n5A"], p["n8F"], p["K58"], p["n58"], p["K85"], p["n85"]
    kR5, kR8, l05, l08, kP5, kP8 = p["kR5"], p["kR8"], p["l05"], p["l08"], p["kP5"], p["kP8"]
    s5, s8  = p["s5"], p["s8"]
    muA, muF = p["muA"], p["muF"]
    du[1] = dP5 = l05 * R5 - kP5 * P5
    du[2] = dP8 = l08 * R8 - kP8 * P8
    du[3] = dR5 = s5*((A/K5A)^n5A)/(1+(P8/K58)^n58+(A/K5A)^n5A) + s05 - kR5 * R5
    du[4] = dR8 = s8*((F/K8F)^n8F)/(1+(F/K8F)^n8F+(P5/K85)^n85) + s08 - kR8 * R8
    du[5] = dA = gA*(maxA * exp(-0.01 * x) - A) + muA
    du[6] = dF = gF*(maxF * exp(-0.01 * (99 - x)) - F) + muF
    du[7] = dr2 = -r2
    du[8] = dc15 = -c15
    du[9] = dc25 = -c25
    du[10] = dc35 = -c35
    du[11] = dr9 = -r9
    du[12] = dc18 = -c18
    du[13] = dc28 = -c28
    du[14] = dc38 = -c38
end

d_ics1 = Dict("R8"=> 0.00900023576123889, "c15"=> 1.761337400332505e-06, "P8"=> 2.9560034347325604e-06, "P5"=> 3.4020485285993286e-07, "c25"=> 9.118941517979407e-05, "c35"=> 0.004721137474792935, "A"=> 1.5e-323, "c18"=> 8.844673354103184e-06, "R5"=> 0.000500023861019239, "c38"=> 0.001867111128611335, "F"=> 1.5e-323, "r2"=> 0.3524443000060298, "r9"=> 0.33001432105772416, "c28"=> 0.0002687664938900425)
svlist = ["P5", "P8", "R5", "R8", "A", "F", "r2", "c15", "c25", "c35", "r9", "c18","c28", "c38"]
u0 = [d_ics1[x] for x in svlist]
tspan = (0.0,1000.0)
p=Dict("x"=>50
    ,"amp"=>0.001
    ,"kR5"=>1.0
    ,"l05"=>1.0
    ,"kP5"=>1.0
    ,"maxA"=>2.3
    ,"maxF"=>2.1
    ,"kR8"=>1.0
    ,"l08"=>0.7
    ,"kP8"=>1.0
    ,"gA"=>1.0
    ,"gF"=>1.0
    ,"gc"=>1.0
    ,"s05"=>0.0
    ,"s5"=>1.0
    ,"s08"=>0.0
    ,"s8"=>1.0
    ,"K5A"=>0.7
    ,"n5A"=>6.0
    ,"K58"=>0.05
    ,"n58"=>6.0
    ,"K85"=>0.18
    ,"n85"=>6.0
    ,"K8F"=>0.3
    ,"n8F"=>6.0
    ,"y"=>1)

