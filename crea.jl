using Nemo, FastArithmetic, BenchmarkTools

const k,u=FiniteField(5,1,"u")
const K,t=PolynomialRing(k,"t")

function create(n::Int64,k::Nemo.Field,c::Int64=5)
  A=Array(fq_nmod,n)
  for i in 1:n
    A[i]=k(rand(0:c-1))
  end
  return A
end

function createPol(n::Int64,K::Nemo.Ring,c::Int64=5)
  k=base_ring(K)
  return K(fq_nmod[k(rand(0:c-1)) for i in 1:n+1])
end

function createIrrPol(n::Int64,K::Nemo.Ring)
  P=createPol(n,K)
  while degree(P)!=n || length(factor(P))!=1
    P=createPol(n,K)
  end
  return P
end

function create2d(m::Int64,n::Int64,k::Nemo.Field)
  A=Array(fq_nmod,(n,m))
  for i in 1:n, j in 1:m
    A[i,j]=k(rand(0:4))
  end
  return A
end
