using ComponentArrays
using ForwardDiff

c_vec = ComponentVector(
    a = [1.0, 2.0, 3.0],
    b = [4.0, 5.0, 6.0]
)

for x in c_vec
    print(typeof(x))
    print(" ")
    println(x)

    x = convert(ForwardDiff.Dual,x)

    print(typeof(x))
    print(" ")
    println(x)

end

y = [convert(ForwardDiff.Dual,x) for x in c_vec]