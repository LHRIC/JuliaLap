using LinearAlgebra

function float_link(initial_float_hdpts, idx1, idx2)
    initial_length = norm(initial_float_hdpts[idx1,:] - initial_float_hdpts[idx2,:])
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        (initial_length - norm(float_hdpts[idx1,:] - float_hdpts[idx2,:])) / initial_length
    end
    return resid_fun
end

function ground_link(initial_float_hdpts, initial_fixed_hdpts, gidx, fidx)
    initial_length = norm(initial_fixed_hdpts[gidx,:] - initial_float_hdpts[fidx,:])
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        (initial_length - norm(fixed_hdpts[gidx] - float_hdpts[fidx])) / initial_length
    end
    return resid_fun
end

function ground_planar(initial_float_hdpts, initial_fixed_hdpts, gidx, fidx, axis)
    resid_fun = (float_hdpts, fixed_hdpts) -> begin
        dot(fixed_hdpts[gidx,:] - float_hdpts[fidx,:], axis)
    return resid_fun
    end 
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

    Rvec = [ground_link(float_hdpts, fixed_hdpts, LIF, LO),     # Lower Fore
            ground_link(float_hdpts, fixed_hdpts, LIA, LO),     # Lower Aft
            ground_link(float_hdpts, fixed_hdpts, UIF, UO),     # Upper Fore
            ground_link(float_hdpts, fixed_hdpts, UIA, UO),     # Upper Aft
            float_link(float_hdpts, OT, IT),                    # Tie Rod
            float_link(float_hdpts, LO, OT),                    # Unsprung - Lower → Tie
            float_link(float_hdpts, UO, OT),                    # Unsprung - Upper → Tie
            float_link(float_hdpts, LO, UO),                    # Unsprung - Lower → Upper
            float_link(float_hdpts, OP, IP),                    # Push/PullRod
            float_link(float_hdpts, IP, SO),                    # Bellcrank - PRod → Shock
            float_link(float_hdpts, IP, BP),                    # Bellcrank - PRod → ARB Pickup
            float_link(float_hdpts, SO, BP),                    # Bellcrank - Shock → ARB Pickup
            ground_link(float_hdpts, fixed_hdpts, BC, BP),      # Bellcrank - Anchor → ARB Pickup
            ground_link(float_hdpts, fixed_hdpts, BC, BP),      # Bellcrank - Anchor → Shock
            float_link(float_hdpts, BP, AP),                    # ARB Droplink
            ground_link(float_hdpts, fixed_hdpts, BC, BP),      # ARB Arm
            ]
    
    bellcrank_axis = cross(float_hdpts[SO,:] - fixed_hdpts[BC,:], float_hdpts[IP,:] - fixed_hdpts[BC,:])
    normalize!(bellcrank_axis)

    push!(Rvec,[
        ground_planar(float_hdpts, fixed_hdpts, BC, IP, bellcrank_axis), # Bellcrank Push/Pull Pickup
        ground_planar(float_hdpts, fixed_hdpts, BC, SO, bellcrank_axis), # Bellcrank Shock Pickup
        ground_planar(float_hdpts, fixed_hdpts, BC, BP, bellcrank_axis), # Bellcrank ARB Pickup
        ground_planar(float_hdpts, fixed_hdpts, BE, AP, [0.0,1.0,0.0])   # ARB Arm
        ]...)

    is_pushrod = norm(float_hdpts[LO,:] - float_hdpts[OP,:]) <= norm(float_hdpts[UO,:] - float_hdpts[OP,:])
    if is_pushrod # Pushrod apex construction
        push!(Rvec,[
            float_link(float_hdpts, LO, OP),                # Apex 1
            ground_link(float_hdpts, fixed_hdpts, LIF, OP), # Apex 2
            ground_link(float_hdpts, fixed_hdpts, LIA, OP)  # Apex 2
        ]...)
    else # Pullrod apex construction
        push!(Rvec,[
            float_link(float_hdpts, UO, OP),                # Apex 1
            ground_link(float_hdpts, fixed_hdpts, UIF, OP), # Apex 2
            ground_link(float_hdpts, fixed_hdpts, UIA, OP)  # Apex 2
        ]...)
    end

    return Rvec
end