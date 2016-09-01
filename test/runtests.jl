#################
#     TESTS
#################

using Nemo, FastArithmetic, Base.Test

function testMonomialDual()
  print("monomialToDual(), dualToMonomial()... ")

  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  P=t^2+t+1
  a0=[k(1),k(0)]
  a1=[k(0),k(1)]
  a2=a0+a1
  c=[k(4),k(2)]

  @test monomialToDual(a0,P)+monomialToDual(a1,P)==monomialToDual(a2,P)

  @test a0==dualToMonomial(monomialToDual(a0,P),P)

  @test a1==dualToMonomial(monomialToDual(a1,P),P)

  @test c==dualToMonomial(monomialToDual(c,P),P)

  @test dualToMonomial(a0,P)+dualToMonomial(a1,P)==dualToMonomial(a2,P)

  @test a0==monomialToDual(dualToMonomial(a0,P),P)

  @test a1==monomialToDual(dualToMonomial(a1,P),P)

  @test c==monomialToDual(dualToMonomial(c,P),P)
  println("PASS")
end

function testMulRem()
  print("mulT(), mulTmid(), remT()... ")

  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  P=t^2-t+2
  Q=2*t-1
  a=[k(2),k(1),k(1),k(0),k(3)]
  b=[k(-1),k(0),k(1),k(2),k(3),k(-1)]

  @test mulT(a,P,2)==t+4

  @test mulTmid(a,P,2)==t+4

  @test mulT(b,Q,4)==4*t^3+3*t^2+2*t+1

  @test mulTmid(b,Q,4)==4*t^3+3*t^2+2*t+1

  l,v=FiniteField(1993,1,"v")
  S,x=PolynomialRing(l,"x")
  fibo=[l(1),l(1)]
  W=x^2-x-1
  Y=t+1

  @test remT(fibo,W,7)==[l(1),l(1),l(2),l(3),l(5),l(8),l(13)]

  @test remT([k(1)],Y,3)==[k(1),k(4),k(1)]

  println("PASS")
end

function testRemTnaif()
  print("remTnaif()... ")
  l,v=FiniteField(1993,1,"v")
  S,x=PolynomialRing(l,"x")
  fibo=[l(1),l(1)]
  W=x^2-x-1

  @test remTnaif(fibo,W,8)==[l(1),l(1),l(2),l(3),l(5),l(8),l(13),l(21)]

  println("PASS")
end

function testBerlekampMassey()
  print("berlekampMassey()... ")

  l,v=FiniteField(1993,1,"v")
  S,x=PolynomialRing(l,"x")

  @test berlekampMassey([l(1),l(1),l(2),l(3),l(5),l(8),l(13),l(21)],2,S)==x^2-x-1

  @test berlekampMassey([l(1),l(1),l(2),l(3)],2,S)==x^2-x-1

  println("PASS")
end

function testEmbedProject()
  print("embed(), project(), computeR()... ")

  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  P=t^2+t+1
  Q=t^3+t+1
  R=computeR(P,Q)
  F=ResidueRing(ZZ,5)
  S,x=PolynomialRing(F,"x")

  @test parent(R)==parent(P)

  @test degree(R)==degree(P)*degree(Q)

  @test length(factor(R))==1

  up=monomialToDual([k(1)],P)
  uq=monomialToDual([k(1)],Q)
  a0=[k(1),k(0)]
  a1=[k(0),k(1)]
  a2=a0+a1
  b0=[k(1),k(0),k(0)]
  b2=[k(0),k(0),k(1)]
  b3=b0+b2

  @test embed(a0,P,uq,Q)+embed(a1,P,uq,Q)==embed(a2,P,uq,Q)

  @test embed(up,P,b0,Q)+embed(up,P,b2,Q)==embed(up,P,b3,Q)

  @test project(dualToMonomial(embed(monomialToDual(a0,P),P,uq,Q),R),P,Q)==a0

  @test project(dualToMonomial(embed(monomialToDual(a1,P),P,uq,Q),R),P,Q)==a1

  @test project(dualToMonomial(embed(up,P,monomialToDual(b0,Q),Q),R),Q,P)==b0

  @test project(dualToMonomial(embed(up,P,monomialToDual(b2,Q),Q),R),Q,P)==b2

  println("PASS")
end

function testPhi1()
  print("phi1(), inversePhi1()... ")

  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  P=t^2+t+1
  Q=t^3+t+1
  b=fq_nmod[k(1) k(2);k(4) k(2); k(1) k(3)]
  a=phi1(b,P,Q)
  R=computeR(P,Q)
  c=dualToMonomial(a,R)
  B=inversePhi1(c,P,Q)

  B=transpose(B)

  for j in 1:3
    B[:,j]=dualToMonomial(B[:,j],P)
  end

  B=transpose(B)

  for j in 1:2
    B[:,j]=monomialToDual(B[:,j],Q)
  end

  @test b==B

  up=monomialToDual([k(1)],P)
  uq=monomialToDual([k(1)],Q)
  r=embed(up,P,uq,Q)
  ur=fq_nmod[k(1) k(0);k(0) k(0); k(0) k(0)]

  for i in 1:2
    ur[:,i]=monomialToDual(ur[:,i],Q)
  end

  d=phi1(ur,P,Q)

  @test r==d

  J=transpose(inversePhi1([k(0),k(1),k(0),k(0),k(0),k(0)],P,Q))

  for i in 1:3
    J[:,i]=dualToMonomial(J[:,i],P)
  end

  @test J==[k(0) k(0) k(0);k(0) k(1) k(0)]

  println("PASS")
end

function testPhi2()
  print("phi2(), inversePhi2()... ")

  k,u=FiniteField(5,1,"u")
  K,t=PolynomialRing(k,"t")
  P=t^2+t+1
  Q=t^3+t+1
  R=computeR(P,Q)
  un=fq_nmod[k(1) k(0) k(0); k(0) k(0) k(0)]
  xy=fq_nmod[k(0) k(0) k(0); k(0) k(1) k(0)]
  b1=fq_nmod[k(2) k(1) k(4); k(3) k(1) k(4)]
  b2=fq_nmod[k(3) k(2) k(3); k(1) k(2) k(0)]

  @test phi2(un,P,Q)==[k(1),k(0),k(0),k(0),k(0),k(0)]

  @test phi2(xy,P,Q)==[k(0),k(1),k(0),k(0),k(0),k(0)]

  @test phi2(b1+b2,P,Q)==phi2(b1,P,Q)+phi2(b2,P,Q)

  b2p=transpose(copy(b2))

  for i in 1:2
    b2p[:,i]=monomialToDual(b2p[:,i],Q)
  end

  @test dualToMonomial(phi1(b2p,P,Q),R)==phi2(b2,P,Q)

  c=inversePhi2(monomialToDual([k(0),k(1),k(0),k(0),k(0),k(0)],R),P,Q)

  for i in 1:3
    c[:,i]=dualToMonomial(c[:,i],P)
  end

  c=transpose(c)

  for i in 1:2
    c[:,i]=dualToMonomial(c[:,i],Q)
  end

  c=transpose(c)

  @test c==xy

  a1=[k(1),k(3),k(2),k(2),k(0),k(1)]
  a2=[k(3),k(3),k(0),k(4),k(1),k(2)]

  @test inversePhi2(a1+a2,P,Q)==inversePhi2(a1,P,Q)+inversePhi2(a2,P,Q)

  q=inversePhi2(a1,P,Q)

  for i in 1:3
    q[:,i]=dualToMonomial(q[:,i],P)
  end

  q=transpose(q)

  for i in 1:2
    q[:,i]=dualToMonomial(q[:,i],Q)
  end

  q=transpose(q)

  @test monomialToDual(phi2(q,P,Q),computeR(P,Q))==a1

  aa=inversePhi2(monomialToDual(phi2(b1,P,Q),R),P,Q)

  for i in 1:3
    aa[:,i]=dualToMonomial(aa[:,i],P)
  end

  aa=transpose(aa)

  for i in 1:2
    aa[:,i]=dualToMonomial(aa[:,i],Q)
  end

  aa=transpose(aa)

  @test aa==b1

  g=inversePhi1(dualToMonomial(a2,R),P,Q)

  for i in 1:2
    g[:,i]=monomialToDual(g[:,i],Q)
  end

  @test transpose(g)==inversePhi2(a2,P,Q)

  println("PASS")
end

function testAll()

  testMonomialDual()
  testMulRem()
  testRemTnaif()
  testBerlekampMassey()
  testEmbedProject()
  testPhi1()
  testPhi2()

  println("\nAll tests passed successfully.\n")

end

testAll()
