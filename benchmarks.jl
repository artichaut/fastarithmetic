include("crea.jl")
include("irreduciblePol.jl")

function benchPhi1()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create2d(j,j+1,k)
    b = @benchmark phi1($a,$P,$Q)
    A[j,1],A[j,2]=j,median(b).time/10^9
  end
  writedlm("phi1.txt",A)
end

function benchinversePhi1()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create(j*(j+1),k)
    b = @benchmark inversePhi1($a,$P,$Q)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("inversePhi1.txt",A)
  end
end

function benchMonomialToDual()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    J = j*(j+1)
    P = createPol(J,K) 
    a = create(J,k)
    b = @benchmark monomialToDual($a,$P)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("monomialToDual.txt",A)
  end
end

function benchDualToMonomial()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    J = j*(j+1)
    P = createPol(J,K)
    a = create(J,k)
    b = @benchmark dualToMonomial($a,$P)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("dualToMonomial.txt",A)
  end
end

function benchEmbed()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    a = create(j,k)
    c = create(j+1,k)
    b = @benchmark embed($a,$P,$c,$Q)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("embed.txt",A)
  end
end

function benchComputeR()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    P = irrPol[j]
    Q = irrPol[j+1]
    b = @benchmark computeR($P,$Q)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("computeR.txt",A)
  end
end

function benchProject()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    a = create(j*(j+1),k)
    P = irrPol[j]
    Q = irrPol[j+1]
    b = @benchmark project($a,$P,$Q)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("project.txt",A)
  end
end

function benchMulTmid()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    a = create(2*j+1,k)
    P = irrPol[j]
    b = @benchmark mulTmid($a,$P,$j)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("mulTmid.txt",A)
  end
end

function benchRemT()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    a = create(j,k)
    P = irrPol[j]
    b = @benchmark remT($a,$P,$j*2)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("remT.txt",A)
  end
end

function benchComputeR()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(20,2))
  for j in 1:20
    println(10*j)
    a = create(10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    b = @benchmark computeR($P,$Q)
    A[j,1],A[j,2]=10*j,median(b).time/10^9
    writedlm("computeR.txt",A)
  end
end

function benchPhi2Inv()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(10,3))
  for j in 1:10
    println(10*j)
    a = create2d(10*j+1,10*j,k)
    abis = create((10*j+1)*10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    R = computeR(P,Q)
    b = @benchmark phi2($a,$P,$Q,$R)
    bbis = @benchmark inversePhi2($abis,$P,$Q,$R)
    A[j,1],A[j,2],A[j,3]=10*j,median(b).time/10^9,median(bbis).time/10^9
    writedlm("phi2inversephi2.txt",A)
  end
end

function benchPhi2()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(20,2))
  for j in 1:20
    println(10*j)
    a = create2d(10*j+1,10*j,k)
    P = irrPol[10*j]
    Q = irrPol[10*j+1]
    R = computeR(P,Q)
    b = @benchmark phi2($a,$P,$Q,$R)
    A[j,1],A[j,2]=10*j,median(b).time/10^9
    writedlm("phi2.txt",A)
  end
end

function benchMulmod()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    P,Q,R=createPol(j*(j+1),K),createPol(j*(j+1),K),createPol(j*(j+1),K)
    b = @benchmark mulmod($P,$Q,$R)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("mulmod.txt",A)
  end
end

function benchMulModT()
  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  A=Array(Float64,(200,2))
  for j in 1:200
    println(j)
    J = j*(j+1)
    a = create(J,k)
    Q,R=createPol(J,K),createPol(J,K)
    b = @benchmark mulModT($a,$Q,$R,$J)
    A[j,1],A[j,2]=j,median(b).time/10^9
    writedlm("mulModT.txt",A)
  end
end
