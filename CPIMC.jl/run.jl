using OnlineStats

include("model.jl")
include("updates.jl")
include("MC.jl")

function main()
  e = Ensemble(200, 2*pi, 0.1, 50.0)
  c = Configuration(Set(),0)

  updates = Set([move_particle])

  measurements =
    [ (Variance(), totalEnergy)
    , (Variance(), particleNumber)
    , (Group([Variance() for i=1:e.cutoff]), occVec)
    ]

  acc = sweep(10^7, 100, updates, measurements, e, c)

  println("accepted/proposed:")
  println("==================")

  for (f,r) in zip(updates,acc)
    println(typeof(f).name.mt.name, "\t(", r[1], "/", r[2], ")")
  end

  println("")
  println("measurements:")
  println("=============")

  for (f,m) in measurements
    if typeof(f) == Variance{Float64,EqualWeight}
      println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f))
    end
  end

  println("occupations:")
  println("============")
  println(mean.(measurements[3][1].stats))
  println("")
  println(std.(measurements[3][1].stats))
  println("")

end

main()
