function hdpt_vec(c_vec)
    # println(typeof(c_vec[:LO]))
    float_hdpts = vcat(
        c_vec[:LO]', # 1 Lower Outboard
        c_vec[:UO]', # 2 Upper Outboard
        c_vec[:OT]', # 3 Outboard Tie
        c_vec[:IT]', # 4 Inboard Tie
        c_vec[:OP]', # 5 Outboard Push/Pull
        c_vec[:IP]', # 6 Inboard Push/Pull
        c_vec[:SO]', # 7 Outboard Shock
        c_vec[:BP]', # 8 Bellcrank ARB Pickup
        c_vec[:AP]'  # 9 ARB Arm
    )

    fixed_hdpts = vcat(
        c_vec[:LIF]', # 1 Lower Inboard Fore
        c_vec[:LIA]', # 2 Lower Inboard Aft
        c_vec[:UIF]', # 3 Upper Inboard Fore
        c_vec[:UIA]', # 4 Upper Inboard Aft
        c_vec[:BC]', # 5 Bellcrank Anchor
        c_vec[:SI]', # 6 Inboard Shock
        c_vec[:BE]' # 7 ARB Bar End
    )

    fixed_hdpts[7,2] = float_hdpts[9,2] # Project ARB axis to arm plane

    return float_hdpts, fixed_hdpts
end