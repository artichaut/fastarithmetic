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
  end
  writedlm("phi1.txt",A)
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
    writedlm("inversePhi1.txt",A)
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
    writedlm("monomialToDual.txt",A)
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
    writedlm("dualToMonomial.txt",A)
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
    writedlm("embed.txt",A)
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
    writedlm("project.txt",A)
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
    writedlm("mulTmid.txt",A)
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
    writedlm("remT.txt",A)
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
    writedlm("computeR.txt",A)
  end
end

function benchPhi2Inv(sizes=1:10)
  A=Array(Float64,(length(sizes),3))
  for (i,j) in enumerate(sizes)
    println(10*j)
    a = create2d(10*j+1,10*j,k)
    abis = create((10*j+1)*10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    R = computeR(P,Q)
    b = @benchmark phi2($a,$P,$Q,$R)
    bbis = @benchmark inversePhi2($abis,$P,$Q,$R)
    A[i,1],A[i,2],A[i,3]=10*j,median(b).time/10^9,median(bbis).time/10^9
    writedlm("phi2inversephi2.txt",A)
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
    writedlm("phi2.txt",A)
  end
end

function benchMulmod(sizes=1:200)
  A=Array(Float64,(length(sizes),2))
  for (i,j) in enumerate(sizes)
    println(j)
    P,Q,R=createPol(j*(j+1),K),createPol(j*(j+1),K),createPol(j*(j+1),K)
    b = @benchmark mulmod($P,$Q,$R)
    A[i,1],A[i,2]=j,median(b).time/10^9
    writedlm("mulmod.txt",A)
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
    writedlm("mulModT.txt",A)
  end
end
