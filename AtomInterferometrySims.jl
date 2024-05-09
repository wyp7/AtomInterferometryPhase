module AtomInterferometrySims

export atom_phase_path_int, Constants, simple_test, test_atom_phase_path_int, init_vel, init_pos

using DifferentialEquations, Parameters, ParameterizedFunctions, LinearAlgebra
using BenchmarkTools, Test, CSV, DataFrames, FastGaussQuadrature, StaticArrays

@with_kw struct Constants
#= Defining necessary physical constants 
These are saved in the struct "constants" as a convenient input for local
functions at the end of the file
=#
    Re_val::Float64 = 6.371e6; # m
    g_val::Float64 = 9.81; # m/s
    theta_lat::Float64 = 42.05; # Evanston; IL
    Omega_val::Float64 = 7.27e-5; # rad/s
    yOmega::Float64 = Omega_val*cosd(theta_lat)
    zOmega::Float64 = Omega_val*sind(theta_lat)
    g = [0, 0, -g_val]
    Re = [0, 0, Re_val]
    Omega = [0, yOmega, zOmega]
    Tzz::Float64 = -2*g_val/Re_val
    Txx::Float64 = g_val/Re_val
    Tyy::Float64 = Txx
    Qzzz::Float64 = 6*g_val/Re_val^2
    Szzzz::Float64 = -24*g_val/Re_val^3
    lambda::Float64 = 689e-9;# m
    k::Float64 = 2*pi/lambda
    hbar::Float64 = 6.62607015e-34/(2*pi); # Js
    m::Float64 = 87.62*1.66*10^-27;#1.45970703e-25; # kg for Sr88
    c::Float64 = 299792458; # m/s
end


function atom_phase_path_int(r0::Vector{Float64}, v0::Vector{Float64}, t0::Float64, T::Float64, n::Real, constants::Constants = Constants())
    
    r0 = SVector{3}(r0) 
    v0 = SVector{3}(v0)

    # Kick velocity
    k_eff = n*constants.k
    v_k = SVector{3,Float64}(0.,0.,constants.hbar*k_eff/constants.m) #350ns

    t, w = gausslegendre(4000)
    # Difference between the accumulated phase between the top arm and the bottom arm from first beamsplitter
    # pulse until the mirror pulse
    SdiffCB, rC, vC, rB, vB = action_diff(r0, r0, v0 + v_k,       v0, SA_F64[  t0, t0 + T], constants,t,w)

    # Phase difference from the mirror pulse until the final beamsplitter pulse
    SdiffED, rE, vE, rD, vD = action_diff(rC, rB, vC - v_k, vB + v_k, SA_F64[t0+T, t0+2*T], constants,t,w)

    # Final beamsplitter pulse
    vE = vE + v_k 

    # Add together to get the propagation phase difference
    phi_prop = SdiffCB + SdiffED

    # Calculate the phase contribution from the laser 
    phi_laser= (phase_laser(r0, k_eff) 
        -       phase_laser(rC, k_eff)
        +       phase_laser(rE, k_eff)
        -       phase_laser(rB, k_eff))

    # Separation phase
    phi_sep = phase_sep(rD, rE, vD, vE, constants.m, constants.Omega, constants.Re, constants.hbar)

    # Add together propagation phase, laser phase, and separation phase for final total phase difference
    phase_output = phi_prop + phi_laser + phi_sep

    return phase_output#, [phi_prop,phi_laser,phi_sep]#rD, rE, phi_sep#pD, pE
end

# use external phase map input
function atom_phase_path_int(r0::Vector{Float64}, v0::Vector{Float64}, t0::Float64, T::Float64, n::Real, phi_map::Array, size::Float64, constants::Constants = Constants())
    #= Calculating atom phase using path integral approach.
    Required inputs: 
    r0: [x, y, z] coordinates
    v0: [vx, vy, vz] coordinates
    t0: start time
    T:  time between pulses
    n:  number of ħk kicks (instantaneous)  
    Optional inputs:
    phi_map: phase map
    size: size of grid for phase map (in meters)
    =#
    
    # Main code
    N = 1023
    #println(1)
    r0 = SA_F64[r0[1],r0[2],r0[3]]
    v0 = SA_F64[v0[1],v0[2],v0[3]]
    # Kick velocity
    k_eff = n*constants.k
    v_k = SA_F64[0., 0., constants.hbar*k_eff/constants.m]

    t, w = gausslegendre(4000)

    # Difference between the accumulated phase between the top arm and the bottom arm from first beamsplitter
    # pulse until the mirror pulse
    SdiffCB, rC, vC, rB, vB = action_diff(r0, r0, v0 + v_k,       v0, SA_F64[  t0, t0 + T], constants,t,w)

    # Phase difference from the mirror pulse until the final beamsplitter pulse
    SdiffED, rE, vE, rD, vD = action_diff(rC, rB, vC - v_k, vB + v_k, SA_F64[t0+T, t0+2*T], constants,t,w)

    # Final beamsplitter pulse
    vE = vE + v_k 

    # Add together to get the propagation phase difference
    phi_prop = SdiffCB + SdiffED

    # Calculate the phase contribution from the laser during the pulses if there is a map
    #=
    phi_laser= (phase_laser(r0,   0, phi_map, size, N, constants, k_eff) 
        -       phase_laser(rC,   T, phi_map, size, N, constants, k_eff)
        +       phase_laser(rE, 2*T, phi_map, size, N, constants, k_eff)
        -       phase_laser(rB,   T, phi_map, size, N, constants, k_eff))
  
    =#
    phi_laser= (phase_laser(r0, k_eff) 
        -       phase_laser(rC, k_eff)
        +       phase_laser(rE, 2*T, phi_map, size, N, constants, k_eff)
        -       phase_laser(rB, k_eff))
    

    # Determine the momentum for the endpoint of each arm
    pD = constants.m * (vD + cross(constants.Omega, rD + constants.Re))
    pE = constants.m * (vE + cross(constants.Omega, rE + constants.Re))

    # Use final momenta to get separation phase
    phi_sep = dot((pD + pE), rD - rE)/(2*constants.hbar)

    # Add together propagation phase, laser phase, and separation phase for final total phase difference
    phase_output = phi_prop + phi_laser + phi_sep

    return phase_output,collect(rE)#, rD, rE, vD, vE
end


# Internal Functions

function phase_sep(rD, rE, vD, vE, m, Omega, Re, hbar) #500ns
    # Determine the momentum for the endpoint of each arm
    pD = m * (vD + cross(Omega, rD + Re))
    pE = m * (vE + cross(Omega, rE + Re))

    # Use final momenta to get separation phase
    phi_sep = 1/(2*hbar) * dot((pD + pE), rD - rE)
    return phi_sep
end


function phase_laser(r::SVector{3,Float64}, k_eff::Float64) #7ns
    # phi_laser assumes pulse is instantaneous. Gravity gradients and finite pulse effects are ignored
    return dot(r,SA_F64[0., 0., k_eff])
end

function phase_laser(r, t, phi_map, size, N, constants, k_eff)
    # phi_laser assumes pulse is instantaneous. Gravity gradients and finite pulse effects are ignored
    # Uses input phase map to get nearest phase value based on atom location
    kvec = [0., 0., k_eff]
    phi0 = get_phi(r, phi_map, size, N)
    out = dot(kvec,r) - k_eff*t/constants.c + phi0
    return out
end

function get_phi(r0, phi_in, size, N)
    
    row_i, col_i, dist_err = find_point(r0, size, N)
    # 'dist_err' shows how far away chosen grid point is from 'r0'

    phi = phi_in[row_i,col_i]
    return phi#, dist_err
end

function find_point(r,size,N)
   # println(1)
    
    # This function finds the closest row & column indices in phase map corresponding to the position given by input 'r'
    grid = range(-size/2,size/2,length=N)
   # println(2)
    dist_err = [0 0]

    col = argmin(abs.(grid.-r[1]))
   # println(3)
    row = argmin(abs.(grid.-r[2]))
    return row, col, dist_err
end


function action_diff(r01::SVector{3,Float64}, r02::SVector{3,Float64}, 
    v01::SVector{3,Float64}, v02::SVector{3,Float64}, 
    tspan::SVector{2,Float64}, constants::Constants,t::Array{Float64},w::Array{Float64})
    #= This function is used for the propagation phase. It outputs the difference instead of the two integrals separately 
    to avoid truncation errors from subtracting two large numbers to get a small difference. =#
    p = SA_F64[constants.yOmega, constants.zOmega, constants.Re_val, constants.g_val, 
    constants.Txx, constants.Tyy, constants.Tzz, constants.Qzzz, constants.Szzzz]
    ff1, _, _, _ , f1 = get_rv(r01,v01,tspan,p,t,w)
    ff2, _, w, c1, f2 = get_rv(r02,v02,tspan,p,t,w)

    action = sum(w.*c1.*eval_LL.(ff1,ff2,Ref(p)))*(constants.m/constants.hbar)
    
    return action, SA_F64[f1[1],f1[2],f1[3]],SA_F64[f1[4],f1[5],f1[6]], SA_F64[f2[1],f2[2],f2[3]],SA_F64[f2[4],f2[5],f2[6]]
end

function eval_LL(f1::SVector{6,Float64},f2::SVector{6,Float64}, p::SVector{9,Float64})
    return L_calc(SA_F64[f1[1],f1[2],f1[3]],SA_F64[f1[4],f1[5],f1[6]],p) - L_calc(SA_F64[f2[1],f2[2],f2[3]],SA_F64[f2[4],f2[5],f2[6]],p)
end

function L_calc(r::SVector{3,Float64},v::SVector{3,Float64},p::SVector{9,Float64})
    vp = v + cross(SA_F64[0., p[1], p[2]],r.+SA_F64[0., 0., p[3]])
    phirRe = (-(dot(SA_F64[0., 0., -p[4]],r) -(1/2)*(p[5]*r[1]^2 + p[6]*r[2]^2 + p[7]*r[3]^2)
    - (1/(3*2))*p[8]*r[3]^3 - 1/(4*3*2)*p[9]*r[3]^4))
    return dot(vp,vp)/2 - phirRe
end

function get_rv(r0::SVector{3,Float64}, v0::SVector{3,Float64}, tspan::SVector{2,Float64},  p::SVector{9,Float64},t::Array{Float64},w::Array{Float64})
    # Sets up and solves ODE to get atom trajectories
    f0 = SVector{6,Float64}([r0;v0])
    
    c1, c2 = (tspan[2]-tspan[1])/2, sum(tspan)/2 #370ns
    t_in = c1.*t.+c2 #8us

    function reduced_EOM(f::SVector{6,Float64}, p::SVector{9,Float64}, t) #18ns
        # Function input for ODE solver in which we reduce the system of three second order ODEs to six first order ODEs
        xpp = (p[1]^2+p[2]^2 - p[5])*f[1] + 2*p[2]*f[5]-2*p[1]*f[6]
    
        ypp = -p[1]*p[2]*p[3] + (p[2]^2 - p[6])*f[2] - p[2]*(p[1]*f[3]+2*f[4])
    
        zpp = (p[1]^2*p[3] - p[1]*p[2]*f[2] - (1/6)*f[3]*(-6*p[1]^2 + 6*p[7] 
                + 3*p[8]*f[3] + p[9]*f[3]^2) - p[4] + 2*p[1]*f[4])
    
        return SA_F64[f[4], f[5], f[6], xpp, ypp, zpp]
    end

    prob = ODEProblem(reduced_EOM,f0,tspan,p,saveat = t_in) #52us
    alg = AutoTsit5(Rosenbrock23())
    function internalnorm(u,t::Float64)
        norm(u)
    end
    sol = solve(prob,alg,internalnorm=internalnorm, abstol=1e-16, reltol=1e-14) #843us
    return sol.u, t_in, w, c1, sol.(tspan[2])
end

function simple_test(tolerance::Float64)

    r0 = [1e-3, 1e-3, 0.]
    v0 = [1e-3, 1e-3, 13.2]

    phase_output, rD, rE, vD, vE = atom_phase_path_int(r0,v0,0.,2.,40)

    true_phase = BigFloat(-1.41023069582194569520700923859e10);
    err = abs((true_phase-phase_output)/true_phase)
    if isapprox(err,0;atol=tolerance)
        printstyled("Simple phase test passed. \n"; bold = true, color =:green)
    else
        printstyled("Simple phase test failed. \n"; bold = true, color =:red)
        println("RelTol = ",err)
    end
end

function rssq(A)
    return sqrt(sum(A.^2))
end

function test_atom_phase_path_int(tolerance::Float64) 
    
    MC_phase = CSV.File("./test/data/MC_phaseout.csv";header = false,types = BigFloat ).Column1
    MC_rD = CSV.File("./test/data/MC_rD.csv";header = false,types = BigFloat ) |> Tables.matrix
    MC_rE = CSV.File("./test/data/MC_rE.csv";header = false,types = BigFloat ) |> Tables.matrix
    MC_vD = CSV.File("./test/data/MC_vD.csv";header = false,types = BigFloat ) |> Tables.matrix
    MC_vE = CSV.File("./test/data/MC_vF.csv";header = false,types = BigFloat ) |> Tables.matrix
    MC_params = CSV.File("./test/data/MC_params.csv";header = false,types = Float64 ) |> Tables.matrix
    MC_names = CSV.File("./test/data/MC_paramnames.csv";header = false,types = String ).Column1

    N = length(MC_phase)

    reltol_phase = Vector{Float64}(undef,N)
    reltol_rD = Vector{Float64}(undef,N)
    reltol_rE = Vector{Float64}(undef,N)
    reltol_vD = Vector{Float64}(undef,N)
    reltol_vE = Vector{Float64}(undef,N)

    for i in 1:N
        r0 = MC_params[i,1:3]
        v0 = MC_params[i,4:6]
        tau = MC_params[i,9]
        n = MC_params[i,10]/Constants().k
        constants = Constants(g_val=MC_params[i,7],Re_val=MC_params[i,8],hbar=MC_params[i,11],m=MC_params[i,12],yOmega=MC_params[i,13],zOmega=MC_params[i,14],
            Txx=MC_params[i,15],Tyy=MC_params[i,16],Tzz=MC_params[i,17],Qzzz=MC_params[i,18],Szzzz=MC_params[i,19])
        phase_output, rD, rE, vD, vE = atom_phase_path_int(r0,v0,0.,tau,n,constants)
        reltol_phase[i] = abs((MC_phase[i]-phase_output)/MC_phase[i])
        reltol_rD[i] = rssq((MC_rD[i,1:3] - rD)./MC_rD[i,1:3])
        reltol_rE[i] = rssq((MC_rE[i,1:3] - rE)./MC_rE[i,1:3])
        reltol_vD[i] = rssq((MC_vD[i,1:3] - vD)./MC_vD[i,1:3])
        reltol_vE[i] = rssq((MC_vE[i,1:3] - vE)./MC_vE[i,1:3])

        if !isapprox(reltol_phase[i],0;atol=tolerance)
            idx = Int64(ceil(i/10))
            if mod(i,10)>0
                idx2 = mod(i,10)
            else
                idx2 = 10
            end
            printstyled(MC_names[idx]," Test Failed at Case ",idx2,"\n"; bold = true, color =:red)
            println("RelTol Phase = ",reltol_phase[i])
            println("RelTol rD = ",reltol_rD[i])
            println("RelTol rE = ",reltol_rE[i])
            println("RelTol vD = ",reltol_vD[i])
            println("RelTol vE = ",reltol_vE[i])
        end
    end

    if isapprox(maximum(reltol_phase),0;atol=tolerance) 
        printstyled("All MC tests passed. \n"; bold = true, color =:green)
    end
    return nothing #reltol_phase
end

function init_vel(T,num)

    # Input must be in K.
    k = 1.38065e-23; # m^2 kg s^-2 K^-1
    m = 1.4431609e-25; # kilograms
    
    # Average velocity determined from temperature [3D]
    v_avg = sqrt(3*k*T/m)
    
    # Linear array of velocities
    v = range(0,v_avg*10,length=num*10)
    
    # probability distribution function
    f = (4*pi*v.^2).*(m/(2*pi*k*T))^(3/2).*(exp.(-m*v.^2/(2*k*T)))
    
    # random 3d velocities normalized to 1
    v_rand = rand(num,3).-0.5
    v_rand = v_rand./norm.(eachrow(v_rand))
    
    # Get probabilities
    P = f./maximum(f)
    
    # Get values of v_mag within probability distribution
    probs = rand(length(f))
    v_mag = v[probs.<P]
    
    # We want to make sure there are enough elements in v_mag for the output.
    if length(v_mag) < num
        probs = rand(length(f))
        v_new = v[probs.<P]
        append!(v_mag, v_new)
    end
    
    # Create 3D velocities
    v_out = v_rand.*v_mag[1:num]
    return v_out, v_mag[1:num]
    
end

function init_pos(R,num)
    # generate random initial positions, centered around [0 0 0]
    # the length of initial position vectors are uniformly distributed in[0 R]
    # where R is the defined radius of the atom cloud
    
    r = 2*(rand(num,3).-0.5)
    r = r./norm.(eachrow(r))
    r_out = r.*rand(num,1).*R
    
    return r_out
end

end #module