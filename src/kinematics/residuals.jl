using LinearAlgebra

function float_link(initial_float_hdpts, idx1, idx2)
    initial_length = norm(initial_float_hdpts[idx1,:] - initial_float_hdpts[idx2,:])
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        p1 = @view float_hdpts[idx1,:]
        p2 = @view float_hdpts[idx2,:]
        # (initial_length - norm(float_hdpts[idx1,:] - float_hdpts[idx2,:])) / initial_length
        (initial_length - norm(p1 - p2)) / initial_length
    end
    return resid_fun
end

function ground_link(initial_float_hdpts, initial_fixed_hdpts, gidx, fidx)
    initial_length = norm(initial_fixed_hdpts[gidx,:] - initial_float_hdpts[fidx,:])
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        p1 = @view fixed_hdpts[gidx,:]
        p2 = @view float_hdpts[fidx,:]
        (initial_length - norm(p1 - p2)) / initial_length
    end
    return resid_fun
end

function ground_planar(gidx, fidx, axis)
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        p1 = @view fixed_hdpts[gidx,:]
        p2 = @view float_hdpts[fidx,:]
        dot(p1 - p2, axis)
    end 
    return resid_fun
end

function float_set(idx, dim, value)
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        p = @view float_hdpts[idx, dim]
        value .- p
    end
    return resid_fun
end

function ctrl_ground_link(gidx, fidx)
    ctrl_fun = (float_hdpts, fixed_hdpts, x) -> begin
        p1 = @view fixed_hdpts[gidx,:]
        p2 = @view  float_hdpts[fidx,:]
        (x - norm(p1 - p2)) / x
    end
    return ctrl_fun
end


function ctrl_float(idx, dim)
    resid_fun = (float_hdpts, fixed_hdpts, x) -> begin
        p = @view float_hdpts[idx, dim]
        x .- p
    end
    return resid_fun
end


function residual_vec(float_hdpts, fixed_hdpts)

    # Floating Hardpoint Indexes
    LO = 1 # Lower Outboard
    UO = 2 # Upper Outboard
    OT = 3 # Outboard Tie
    IT = 4 # Inboard Tie
    OP = 5 # Outboard Push/Pull
    IP = 6 # Inboard Push/Pull
    SO = 7 # Shock Outboard
    BP = 8 # ARB Bellcrank Pickup
    AP = 9 # ARB Arm Pickup

    # Grounded Hardpoint Indexes
    LIF = 1 # Lower Inboard Fore
    LIA = 2 # Lower Inboard Aft
    UIF = 3 # Upper Inboard Fore
    UIA = 4 # Upper Inboard Aft
    BC = 5  # Bellcrank Anchor
    SI = 6  # Shock Inboard
    BE = 7  # ARB Bar End

    Rvec = [
            ground_link(float_hdpts, fixed_hdpts, LIF, LO),     # 1 Lower Fore
            ground_link(float_hdpts, fixed_hdpts, LIA, LO),     # 2 Lower Aft
            ground_link(float_hdpts, fixed_hdpts, UIF, UO),     # 3 Upper Fore
            ground_link(float_hdpts, fixed_hdpts, UIA, UO),     # 4 Upper Aft
            float_link(float_hdpts, OT, IT),                    # 5 Tie Rod
            float_link(float_hdpts, LO, OT),                    # 6 Unsprung - Lower → Tie
            float_link(float_hdpts, UO, OT),                    # 7 Unsprung - Upper → Tie
            float_link(float_hdpts, LO, UO),                    # 8 Unsprung - Lower → Upper
            float_link(float_hdpts, OP, IP),                    # 9 Push/PullRod
            float_link(float_hdpts, IP, SO),                    # 10 Bellcrank - PRod → Shock
            float_link(float_hdpts, IP, BP),                    # 11 Bellcrank - PRod → ARB Pickup
            float_link(float_hdpts, SO, BP),                    # 12 Bellcrank - Shock → ARB Pickup
            ground_link(float_hdpts, fixed_hdpts, BC, BP),      # 13 Bellcrank - Anchor → ARB Pickup
            ground_link(float_hdpts, fixed_hdpts, BC, SO),      # 14 Bellcrank - Anchor → Shock
            float_link(float_hdpts, BP, AP),                    # 15 ARB Droplink
            ground_link(float_hdpts, fixed_hdpts, BE, AP)       # 16 ARB Arm
            ]
    
    bellcrank_axis = cross(float_hdpts[SO,:] - fixed_hdpts[BC,:], float_hdpts[IP,:] - fixed_hdpts[BC,:])
    normalize!(bellcrank_axis)

    push!(Rvec,[
        ground_planar(BC, IP, bellcrank_axis), # 17 Bellcrank Push/Pull Pickup
        ground_planar(BC, SO, bellcrank_axis), # 18 Bellcrank Shock Pickup
        ground_planar(BC, BP, bellcrank_axis), # 19 Bellcrank ARB Pickup
        ground_planar(BE, AP, [0.0,1.0,0.0])   # 20 ARB Arm
        ]...)

    is_pushrod = norm(float_hdpts[LO,:] - float_hdpts[OP,:]) <= norm(float_hdpts[UO,:] - float_hdpts[OP,:])
    if is_pushrod # Pushrod apex construction
        push!(Rvec,[
            float_link(float_hdpts, LO, OP),                # 21 Apex 1
            ground_link(float_hdpts, fixed_hdpts, LIF, OP), # 22 Apex 2
            ground_link(float_hdpts, fixed_hdpts, LIA, OP)  # 23 Apex 2
        ]...)
    else # Pullrod apex construction
        push!(Rvec,[
            float_link(float_hdpts, UO, OP),                # 21 Apex 1
            ground_link(float_hdpts, fixed_hdpts, UIF, OP), # 22 Apex 2
            ground_link(float_hdpts, fixed_hdpts, UIA, OP)  # 23 Apex 2
        ]...)
    end

    push!(Rvec,[
            float_set(IT, 1, float_hdpts[IT,1]),    # 24 Set Outboard Tie X
            float_set(IT, 3, float_hdpts[IT,3])     # 25 Set Outboard Tie Z
        ]...)

    Cvec = [
            ctrl_ground_link(SI, SO),   # 26 Control Shock Length
            ctrl_float(IT, 2)           # 27 Control Outboard Tie Y-pos
            ]

    return Rvec, Cvec
end