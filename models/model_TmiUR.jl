function model(du,u,p,t)
    ts = [50.0, 100.0, 200.0, 900.0]
    freq = 0.2
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
    P5, P8, R5, R8, A, F, r2, c15, c25, c35, r9, c18, c28, c38, c48 = u
    x, gA, gF, maxA, maxF, s05, s08, K5A, K8F, n5A, n8F, K58, n58  = p["x"], p["gA"], p["gF"], p["maxA"], p["maxF"], p["s05"], p["s08"], p["K5A"], p["K8F"], p["n5A"], p["n8F"], p["K58"], p["n58"]
    K9F, K2A, n9F, n2A, K25, K98, n25, n98 = p["K9F"], p["K2A"], p["n9F"], p["n2A"], p["K25"], p["K98"], p["n25"], p["n98"]
    gc, K = p["gc"], p["K"]
    l05, l15, l25, l35, kP5 = p["l05"], p["l15"], p["l25"], p["l35"], p["kP5"]
    l08, l18, l28, l38, l48, kP8 = p["l08"], p["l18"], p["l28"], p["l38"], p["l48"], p["kP8"]
    kR5, kR8, s2, s9 = p["kR5"], p["kR8"], p["s2"], p["s9"]
    kr2, kr9  = p["kr2"], p["kr9"]
    s5, s8  = p["s5"], p["s8"]
    muA, muF = p["muA"], p["muF"]
    b15 = b25 = b35 = p["b5"]
    a15 = a25 = a35 = p["a5"]
    b18 = b28 = b38 = b48 = p["b8"]
    a18 = a28 = a38 = a48 = p["a8"]
    du[1] = dP5 = l05 * (R5 - 3 * c15 - 3 * c25 - c35) + l15 * 3 * c15 + l25 * 3 * c25 + l35 * c35 - kP5 * P5
    du[2] = dP8 = l08 * (R8 - 4 * c18 - 6 * c28 - 4 * c38 - c48) + l18 * 4 * c18 + l28 * 6 * c28 + 4 * l38 * c38 + l48 * c48 - kP8 * P8
    du[3] = dR5  = s5*((A/K5A)^n5A)/(1+(P8/K58)^n58+(A/K5A)^n5A) + s05 - kR5 * ((R5- 3 * c15 - 3 * c25 - 1 * c35 ) + a15 * 3 * c15 + a25 * 3 * c25 + a35 * 1 * c35)
    du[4] = dR8 = s8*((F/K8F)^n8F)/(1+(F/K8F)^n8F) + s08 - kR8 * ((R8 - 4 * c18 - 6 * c28 - 4 * c38 - c48) + a18 * 4 * c18 + a28 * 6 * c28 + a38 * 4 * c38 + a48 * c48)
    du[5] = dA = gA*(maxA * exp(-0.01 * x) - A) + muA
    du[6] = dF = gF*(maxF * exp(-0.01 * (99 - x)) - F) + muF
    du[7] = dr2 = s2 / (1 + (A/K2A)^n2A + (P5/K25)^n25) - kr2 * (r2 - 1 * 3 * c15 - 2 * 3 * c25 - 3 * 1 * c35 ) - b15 * 1 * 3 * c15 - b25 * 2 * 3 * c25 - b35 * 3 * 1 * c35
    du[8] = dc15 = gc*(K*(R5-3*c15-3*c25-c35)*(r2-3*c15-3*2*c25-3*c35)-c15)
    du[9] = dc25 = gc * (K * c15 * (r2 - 3 * c15 - 3 * 2 * c25 - 3 * c35) - c25)
    du[10] = dc35 = gc * (K * c25 * (r2 - 3 * c15 - 3 * 2 * c25 - 3 * c35) - c35)
    du[11] = dr9 = s9 - kr9 * (r9 - 1 * 4 * c18 - 2 * 6 * c28 - 3 * 4 * c38 - 4 * 1 * c48) - b18 * 1 * 4 * c18  - b28 * 2 * 6 * c28 - b38 * 3 * 4 * c38 - b48 * 4 * 1 * c48
    du[12] = dc18 = gc * (K * (R8 - 4 * c18 - 6 * c28 - 4 * c38 - c48) * (r9 - 1 * 4 * c18 - 2 * 6 * c28 - 3 * 4 * c38 - 4 * 1 * c48) - c18)
    du[13] = dc28 = gc * (K * c18 * (r9 - 1 * 4 * c18 - 2 * 6 * c28 - 3 * 4 * c38 - 4 * 1 * c48) - c28)
    du[14] = dc38 = gc * (K * c28 * (r9 - 1 * 4 * c18 - 2 * 6 * c28 - 3 * 4 * c38 - 4 * 1 * c48) - c38)
    du[15] = dc48 = gc * (K * c38 * (r9 - 1 * 4 * c18 - 2 * 6 * c28 - 3 * 4 * c38 - 4 * 1 * c48) - c48)
end


d_ics1 = Dict("R8"=> 0.00900023576123889, "c15"=> 1.761337400332505e-06, "P8"=> 2.9560034347325604e-06, "P5"=> 3.4020485285993286e-07, "c25"=> 9.118941517979407e-05, "c35"=> 0.004721137474792935, "A"=> 1.5e-323, "c18"=> 8.844673354103184e-06, "R5"=> 0.000500023861019239, "c38"=> 0.001867111128611335, "F"=> 1.5e-323, "r2"=> 0.3524443000060298, "r9"=> 0.33001432105772416, "c28"=> 0.0002687664938900425, "c48"=>0.001)
svlist = ["P5", "P8", "R5", "R8", "A", "F", "r2", "c15", "c25", "c35", "r9", "c18","c28", "c38", "c48"]
u0 = [d_ics1[x] for x in svlist]
tspan = (0.0,1000.0)
p=Dict("x"=>50
    ,"amp"=>0.001
    ,"s2"=>1.0
    ,"kR5"=>1.0
    ,"kr2"=>1.0
    ,"b5"=>1.0
    ,"a5"=>1.0
    ,"b8"=>1.0
    ,"a8"=>1.0
    ,"l05"=>1.0
    ,"l15"=>0.010
    ,"l25"=>0.010
    ,"l35"=>0.010
    ,"kP5"=>1.0
    ,"K"=>1000
    ,"s9"=>0.5
    ,"kR8"=>1.0
    ,"kR8"=>1.0
    ,"kr9"=>1.0
    ,"l08"=>1.0
    ,"l18"=>0.010
    ,"l28"=>0.010
    ,"l38"=>0.010
    ,"l48"=>0.010
    ,"kP8"=>1.0
    ,"gA"=>1.0
    ,"gF"=>1.0
    ,"gc"=>10.0
    ,"s05"=>0.15
    ,"s5"=>1.0
    ,"s08"=>0.05
    ,"s8"=>1.0
    ,"K5A"=>0.98
    ,"n5A"=>6.0
    ,"K58"=>0.1
    ,"n58"=>40.0
    ,"K8F"=>0.95
    ,"n8F"=>40.0
    ,"K9F"=>1000
    ,"n9F"=>6.0
    ,"K2A"=>0.90
    ,"n2A"=>6.0
    ,"K25"=>100
    ,"n25"=>8.0
    ,"K98"=>1000
    ,"n98"=>8.0
    ,"y"=>1)
