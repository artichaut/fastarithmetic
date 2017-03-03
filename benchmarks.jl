include("crea.jl")
include("irreduciblePol.jl")

function benchPhi1(sizes=2:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)      # fails for fmpz if j=1
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create2d(j,j+1,k)
    b = @benchmark phi1($a,$P,$Q)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/phi1.txt",A)
  end
end

function benchPhi1_pre(sizes=2:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)      # fails for fmpz if j=1
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create2d(j,j+1,k)
    up = monomialToDual([k(1)],P)
    b = @benchmark phi1_pre($a,$P,$Q, $up)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/phi1_pre.txt",A)
  end
end

function benchinversePhi1(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create(j*(j+1),k)
    b = @benchmark inversePhi1($a,$P,$Q)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/inversePhi1.txt",A)
  end
end

function benchinversePhi1_pre(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create(j*(j+1),k)
    up = monomialToDual(T[k(1)],P)
    b = @benchmark inversePhi1_pre($a,$P,$Q,$up)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/inversePhi1_pre.txt",A)
  end
end

function benchMonomialToDual(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    J = j*(j+1)
    P = createPol(J,K)
    a = create(J,k)
    b = @benchmark monomialToDual($a,$P)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/monomialToDual.txt",A)
  end
end

function benchMonomialToDual_pre(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    J = j*(j+1)
    P = createPol(J,K)
    m = degree(P)
    TP = reverse(P, m+1)
    TP = gcdinv(TP, t^m)[2] 
    a = create(J,k)
    b = @benchmark monomialToDual_pre($a,$P, $(TP))
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/monomialToDual_pre.txt",A)
  end
end

function benchDualToMonomial(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    J = j*(j+1)
    P = createPol(J,K)
    a = create(J,k)
    b = @benchmark dualToMonomial($a,$P)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/dualToMonomial.txt",A)
  end
end

function benchDualToMonomial_pre(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    J = j*(j+1)
    P = createPol(J,K)
    S = gcdinv(derivative(P),P)[2]
    a = create(J,k)
    b = @benchmark dualToMonomial_pre($a,$P,$S)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/dualToMonomial_pre.txt",A)
  end
end

function benchEmbed(sizes=2:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)      # fails for fmpz if j=1
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create(j,k)
    c = create(j+1,k)
    b = @benchmark embed($a,$P,$c,$Q)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/embed.txt",A)
  end
end

function benchProject(sizes=2:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)      # fails for fmpz if j=1
    println(j)
    a = create(j*(j+1),k)
    P = irrPol[j]
    Q = irrPol[j+1]
    b = @benchmark project($a,$P,$Q)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/project.txt",A)
  end
end

function benchMulTmid(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    a = create(2*j+1,k)
    P = irrPol[j]
    b = @benchmark mulTmid($a,$P,$j)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/mulTmid.txt",A)
  end
end

function benchRemT(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    a = create(j,k)
    P = irrPol[j]
    b = @benchmark remT($a,$P,$j*2)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/remT.txt",A)
  end
end

function benchComputeR(sizes=1:20)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(10*j)
    a = create(10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    b = @benchmark computeR($P,$Q)
    A[i,1],A[i,2]=10*j,median(b).time/10^9
    writedlm("benchmarks/computeR.txt",A)
  end
end

function benchPhi2InvCompR(sizes=2:200)
  A=Array(Float64,(length(sizes),4))
  for (i,j) in enumerate(sizes)
    println(j)
    a = create2d(j+1,j,k)
    abis = create(j*(j+1),k)
    P = irrPol[j]
    Q = irrPol[j+1]
    up = monomialToDual(T[k(1)],P)
    R = @timed computeR(P,Q)
    b = @benchmark phi2_pre($a,$P,$Q,$(R[1]), $up)
    bbis = @benchmark inversePhi2($abis,$P,$Q,$(R[1]), $up)
    A[i,1],A[i,2],A[i,3],A[i,4]=j,median(b).time/10^9,median(bbis).time/10^9,R[2]
    writedlm("benchmarks/phi2inversephi2CompR_pre.txt",A)
  end
end

function benchPhi2InvCompR_pre(sizes=2:200)
  A=Array(Float64,(length(sizes),4))
  for (i,j) in enumerate(sizes)
    println(j)
    a = create2d(j+1,j,k)
    abis = create(j*(j+1),k)
    P = irrPol[j]
    Q = irrPol[j+1]
    up = monomialToDual(T[k(1)],P)
    R = @timed computeR(P,Q)
    b = @benchmark phi2_pre($a,$P,$Q,$(R[1]), $up)
    bbis = @benchmark inversePhi2_pre($abis,$P,$Q,$(R[1]), $up)
    A[i,1],A[i,2],A[i,3],A[i,4]=j,median(b).time/10^9,median(bbis).time/10^9,R[2]
    writedlm("benchmarks/phi2inversephi2CompR_pre.txt",A)
  end
end

function benchPhi2(sizes=1:20)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(10*j)
    a = create2d(10*j+1,10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    R = computeR(P,Q)
    b = @benchmark phi2($a,$P,$Q,$R)
    A[i,1],A[i,2]=10*j,median(b).time/10^9
    writedlm("benchmarks/phi2.txt",A)
  end
end

function benchMulmod(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    P,Q,R=createPol(j*(j+1),K),createPol(j*(j+1),K),createPol(j*(j+1),K)
    b = @benchmark mulmod($P,$Q,$R)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/mulmod.txt",A)
  end
end

function benchMulModT(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    J = j*(j+1)
    a = create(J,k)
    Q,R=createPol(J,K),createPol(J,K)
    b = @benchmark mulModT($a,$Q,$R,$J)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("benchmarks/mulModT.txt",A)
  end
end
