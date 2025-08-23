using ComponentArrays

function kin_params(c_array)
        
        float_hdpts = c_array[(
            :LO, # Lower Outboard
            :UO, # Upper Outboard
            :OT, # Outboard Tie
            :IT,
            :OP, # Outboard Push/Pull
            :IP, # Inboard Push/Pull
            :SO, # Outboard Shock
            :BP, # Bellcrank ARB Pickup
            :AP  # ARB Arm
        )]

        fixed_hdpts = c_array[(
            :LIF, # Lower Inboard Fore
            :LIA, # Lower Inboard Aft
            :UIF, # Upper Inboard Fore
            :UIA, # Upper Inboard Aft
            :BC,  # Bellcrank Anchor
            :SI,  # Inboard Shock
            :BE   # ARB Bar End
        )]

        fixed_hdpts.BE[2] = float_hdpts[:AP][2]

        # ctrl_hdpts = c_array[(
        #     :IT, # Inboard Tie Rod
        # )]

        hdpts = kinHdpts(
            hdptWrapper(float_hdpts), 
            hdptWrapper(fixed_hdpts))
            # hdptWrapper(ctrl_hdpts))

        residual_vector = Vector{Union{Linkage, Revolute, Linear, Forced}}([
            Linkage(hdpts.float_hdpts, hdpts.fixed_hdpts, :LO, :LIF),   # Lower Fore A-Arm
            Linkage(hdpts.float_hdpts, hdpts.fixed_hdpts, :LO, :LIA),   # Lower Aft A-Arm
            Linkage(hdpts.float_hdpts, hdpts.fixed_hdpts, :UO, :UIF),   # Upper Fore A-Arm
            Linkage(hdpts.float_hdpts, hdpts.fixed_hdpts, :UO, :UIA),   # Upper Aft A-Arm
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :OT, :IT),    # Tie Rod
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :OT, :LO),    # Unsprung - Tie → Lower
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :OT, :UO),    # Unsprung - Tie → Upper
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :LO, :UO),    # Unsprung - Lower → Upper
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :OP, :IP),    # Push/Pull Rod
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :IP, :SO),    # Bellcrank - PRod → Shock 
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :IP, :BP),    # Bellcrank - PRod → ARB Pickup
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :BP, :SO),    # Bellcrank - ARB Pickup → Shock
            Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :BC, :BP),    # Bellcrank - Anchor → ARB Pickup
            Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :BC, :SO),    # Bellcrank - Anchor → Shock
            Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :BP, :AP),    # ARB Droplink
            Linkage(hdpts.float_hdpts, hdpts.fixed_hdpts, :AP, :BE),    # ARB Arm
        ])

        bellcrank_axis = cross(
            hdpts.float_hdpts[:SO] - hdpts.fixed_hdpts[:BC],
            hdpts.float_hdpts[:IP] - hdpts.fixed_hdpts[:BC]
        )
        normalize!(bellcrank_axis)
        
        push!(residual_vector,[
            Revolute(hdpts.fixed_hdpts, hdpts.float_hdpts, :BC, :IP, bellcrank_axis),   # Bellcrank Plane - PRod
            Revolute(hdpts.fixed_hdpts, hdpts.float_hdpts, :BC, :SO, bellcrank_axis),   # Bellcrank Plane - Shock
            Revolute(hdpts.fixed_hdpts, hdpts.float_hdpts, :BC, :BP, bellcrank_axis),   # Bellcrank Plane - ARB Pickup
        ]...)

        push!(residual_vector,[
            Revolute(hdpts.fixed_hdpts, hdpts.float_hdpts, :BE, :AP, [0.0, 1.0, 0.0])   # ARB Plane #17
        ]...)

        
        is_pushrod = norm(hdpts.float_hdpts[:LO] - hdpts.float_hdpts[:OP]) <= norm(hdpts.float_hdpts[:UO] - hdpts.float_hdpts[:OP])
        if is_pushrod
            # println("PUSHROD!!!")
            push!(residual_vector,[
                Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :LO, :OP), # Apex - Lower Outboard → PRod #18
                Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :LIF, :OP), # Apex - Lower Inboard Fore → PRod #19
                Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :LIA, :OP), # Apex - Lower Inboard Aft → PRod #20
                ]...)
        else
            # println("PULLROD!!!")
            push!(residual_vector,[
                Linkage(hdpts.float_hdpts, hdpts.float_hdpts, :UO, :OP), # Apex - Upper Outboard → PRod #21
                Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :UIF, :OP), # Apex - Upper Inboard Fore → PRod #22
                Linkage(hdpts.fixed_hdpts, hdpts.float_hdpts, :UIA, :OP), # Apex - Upper Inboard Aft → PRod #23
                ]...)
        end
                
        initial_shock_length = norm(hdpts.float_hdpts[:SO] - hdpts.fixed_hdpts[:SI])
                
        shock_ctrl = Linear(hdpts.float_hdpts, hdpts.fixed_hdpts, :SO, :SI, initial_shock_length)
        push!(residual_vector,[
            shock_ctrl
        ]...)

        steer_ctrl = Forced(hdpts.float_hdpts, :IT, hdpts.float_hdpts[:IT][2], 2)
        push!(residual_vector,[
            steer_ctrl,                                                     # Controllable Outboard Tie Y
            Forced(hdpts.float_hdpts, :IT, hdpts.float_hdpts[:IT][1], 1),   # Fixed Outboard Tie X
            Forced(hdpts.float_hdpts, :IT, hdpts.float_hdpts[:IT][3], 3),   # Fixed Outboard Tie Z
        ]...)
        
    return kinParameters(hdpts, residual_vector, shock_ctrl, steer_ctrl)
end

