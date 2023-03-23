using PBWDeformations
using Oscar
using BenchmarkTools


printdivider() = println('-'^25)
printdividerdouble() = println('='^25)

println("A1, [1]")
sp, _ = smash_product_lie(QQ, 'A', 1, [1])
for deg in 0:10
    @btime all_pbwdeformations(sp, $deg)
end
printdivider()
println("A3, [1,0,0]")
sp, _ = smash_product_lie(QQ, 'A', 3, [1, 0, 0])
for deg in 0:2
    @btime all_pbwdeformations(sp, $deg)
end
printdividerdouble()

println("B2, [1,0]")
sp, _ = smash_product_lie(QQ, 'B', 2, [1, 0])
for deg in 0:2
    @btime all_pbwdeformations(sp, $deg)
end
printdividerdouble()

println("B2, [0,1]")
sp, _ = smash_product_lie(QQ, 'C', 2, [0, 1])
for deg in 0:1
    @btime all_pbwdeformations(sp, $deg)
end
printdividerdouble()

println("D4, [1,0,0,0]")
sp, _ = smash_product_lie(QQ, 'D', 4, [1, 0, 0, 0])
for deg in 0:1
    @btime all_pbwdeformations(sp, $deg)
end
printdivider()
println("D4, [0,1,0,0]")
sp, _ = smash_product_lie(QQ, 'D', 4, [0, 1, 0, 0])
for deg in 0:0
    @btime all_pbwdeformations(sp, $deg)
end
printdivider()
println("D4, [0,0,1,0]")
sp, _ = smash_product_lie(QQ, 'D', 4, [0, 0, 1, 0])
for deg in 0:0
    @btime all_pbwdeformations(sp, $deg)
end
printdivider()
println("D4, [0,0,0,1]")
sp, _ = smash_product_lie(QQ, 'D', 4, [0, 0, 0, 1])
for deg in 0:1
    @btime all_pbwdeformations(sp, $deg)
end
