using Nemo, FastArithmetic, BenchmarkTools

const k=ResidueRing(ZZ,5) # ,u=FiniteField(5,1,"u")
const K,t=PolynomialRing(k,"t")

Base.size(k::Nemo.GenResRing{Nemo.fmpz}) = Int(modulus(k))
Base.size{T<:Nemo.FinField}(k::T) = Int(order(k))

function create(n::Int64,k::Nemo.Ring)
  A=Array(k,n)
  c=size(k)-1
  for i in 1:n
    A[i]=k(rand(0:c))
  end
  return A
end

function createPol(n::Int64,K::Nemo.Ring)
  k=base_ring(K)
  return K(create(n,k))
end

function createIrrPol(n::Int64,K::Nemo.Ring)
  P=createPol(n,K)
  while degree(P)!=n || length(factor(P))!=1
    P=createPol(n,K)
  end
  return P
end

function create2d(m::Int64,n::Int64,k::Nemo.Ring)
  A=Array(k,n,m)
  c=size(k)-1
  for i in 1:n, j in 1:m
    A[i,j]=k(rand(0:c))
  end
  return A
end
