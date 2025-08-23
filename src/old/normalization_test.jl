using StaticArrays
using LinearAlgebra
using BenchmarkTools

vector1 = SVector(1.0,2.0,3.0,4.0)
vector2 = MVector(1.0,2.0,3.0,4.0)
vector3 = [1.0,2.0,3.0,4.0]

matrix1 = @MMatrix [1.0 2.0 3.0 ; 4.0 5.0 6.0 ; 7.0 8.0 9.0]
matrix2 = @SMatrix [1.0 2.0 3.0 ; 4.0 5.0 6.0 ; 7.0 8.0 9.0]
matrix3 = [1.0 2.0 3.0 ; 4.0 5.0 6.0 ; 7.0 8.0 9.0]

@btime A = $transpose(matrix1)
@btime B = $transpose(matrix2)

print(matrix1)

nothing