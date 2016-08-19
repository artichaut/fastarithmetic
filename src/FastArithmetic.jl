module FastArithmetic

using Nemo

export monomialToDual

"""
    monomialToDual{T <: FieldElem}(a::Array{T,1},P::PolyElem{T})

Convert monomial (canonical) coordinates "a" to dual coordinates
with respect to the trace generated by P.

# Remark
* a and P are over a field k which must be a perfect field ;
* P must be monic and squarefree ;
* k[x]/(P) is not necessarily a field since P is not necessarily irreducible.
"""
function monomialToDual{T <: FieldElem}(a::Array{T,1},P::PolyElem{T})
  k::Nemo.Field=parent(a[1])
  @assert k==base_ring(P)
  m::Int64=degree(P)
  R::Nemo.Ring=parent(P)
  t::PolyElem{T}=gen(R)
  Q::PolyElem{T}=reverse(P,m+1)
  Q=gcdinv(Q,t^m)[2]
  A::PolyElem{T}=R(a)
  b::PolyElem{T}=(reverse((derivative(P)*A)%P,m)*Q)%(t^m)
  return T[coeff(b,i) for i in 0:(m-1)]
end

export dualToMonomial

"""
    dualToMonomial{T}(b::Array{T,1},P::PolyElem{T})

Convert dual coordinates "b" to monomial coordinates with respect to the trace
generated by P.

# Remark
* a and P are over a field k which must be a perfect field ;
* P must be monic and squarefree ;
* k[x]/(P) is not necessarily a field since P is not necessarily irreducible.
"""
function dualToMonomial{T}(b::Array{T,1},P::PolyElem{T})
  k::Nemo.Field=parent(b[1])
  @assert k==base_ring(P)
  m::Int64=degree(P)
  R::Nemo.Ring=parent(P)
  t::PolyElem{T}=gen(R)
  S::PolyElem{T}=gcdinv(derivative(P),P)[2]
  c::PolyElem{T}=(reverse(P,m+1)*R(b))%(t^m)
  c=reverse(c,m)
  d::PolyElem{T}=(c*S)%P
  return T[coeff(d,i) for i in 0:(m-1)]
end

export mulT

"""
    mulT{T}(c::Array{T,1},P::PolyElem{T},n::Int64)

The tranposition of the algorithm of multiplication by P.

# Arguments
* c::Array{T,1} must have length m (degree of P) + n.

# Remark
* closely linked with the middle product [1] ;
* I don't really see why n is an argument, since we could obtain it by
computing n = length(c) - degree(P).

# References
* [1] : G. Hanrot, M. Quercia, and P. Zimmerman. The middle product algorithm I.
Appl. Algebra Engrg. Comm. Comput., 14(6):415-438, 2004.
"""
function mulT{T}(c::Array{T,1},P::PolyElem{T},n::Int64)
  m::Int64=degree(P)
  k::Nemo.Field=base_ring(P)
  @assert k==parent(c[1])
  R::Nemo.Ring=parent(P)
  t::PolyElem{T}=gen(R)
  p::Array{T,1}=fq_nmod[k(coeff(P,j)) for j in 0:m]
  b::Array{T,1}=Array(fq_nmod,n+1) # array filled with "#undef"
  for i in 1:(n+1)
    b[i]=k(0)
  end
  for i in range(m+n,-1,m+n+1)
    for j in range(min(m,i),-1,min(m,i)-max(0,i-n)+1)
      b[i-j+1]=b[i-j+1]+p[j+1]*c[i+1]
    end
  end
  return R(fq_nmod[b[i] for i in 1:(n+1)])
end

export mulTmid

"""
    mulTmid{T}(c::Array{T,1},P::PolyElem{T},n::Int64)

The tranposition of the algorithm of multiplication by P. Using middle product.

# Arguments
* c::Array{T,1} must have length m (degree of P) + n.

# Remark
* the middle product formula is in [1]
* I don't really see why n is an argument, since we could obtain it by
computing n = length(c) - degree(P).

# References
* [1] : G. Hanrot, M. Quercia, and P. Zimmerman. The middle product algorithm I.
Appl. Algebra Engrg. Comm. Comput., 14(6):415-438, 2004.
"""
function mulTmid{T}(c::Array{T,1},P::PolyElem{T},n::Int64)
  m::Int64=degree(P)
  k::Nemo.Field=base_ring(P)
  @assert k==parent(c[1])
  R::Nemo.Ring=parent(P)
  t::PolyElem{T}=gen(R)
  C::PolyElem{T}=R(c)
  Q::PolyElem{T}=reverse(P,m+1)
  prod::PolyElem{T}=Q*C
  prod=prod%t^(n+m+1)
  prod=R(T[coeff(prod,j+m) for j in 0:degree(prod)]) # we're looking a bit too factor
  return prod
end

export remT

"""
    remT{T}(r::Array{T,1},P::PolyElem{T},n::Int64)

Transposition of the remainder by P algorithm.

An other linear extension algorithm. Take the r first values of a linear
recurring sequence of minimal polynomial P and compute the n+1 first values.

# Output
A polynomial of degree n which coefficients are the values of the sequence.
"""
function remT{T}(r::Array{T,1},P::PolyElem{T},n::Int64)
  m::Int64=degree(P)
  k::Nemo.Field=base_ring(P)
  @assert k==parent(r[1])
  K::Nemo.Ring=parent(P)
  t::PolyElem{T}=gen(K)
  R::PolyElem{T}=K(fq_nmod[r[i] for i in 1:m])# useless creation of a list... memory !
  α::PolyElem{T}=reverse(P,m+1)
  α=gcdinv(α,t^(n-m+1))[2]
  b::Array{T,1}=copy(r)
  while length(b)<(n+1)
    push!(b,k(0))
  end
  return R-(t^m)*((α*mulT(b,P,n-m))%(t^(n-m+1)))
end

export remTnaif

"""
    remTnaif{T}(r::Array{T,1},P::PolyElem{T},n::Int64)

Linear extension algorithm.

Take the m first elements of a Linear recurring sequence with minimal
polynomial P and compute the n first.
"""
function remTnaif{T}(r::Array{T,1},P::PolyElem{T},n::Int64)
  m::Int64=degree(P)
  p::Array{T,1}=T[coeff(P,j) for j in 0:m]
  b::Array{T,1}=copy(r)
  while length(b)<n
      s=sum([-1*p[j]*b[end-m+j] for j in 1:m])
      push!(b,p[end]^(-1)*s)
  end
  return b
end

export embed

"""
    embed{T}(b::Array{T,1},P::PolyElem{T},c::Array{T,1},Q::PolyElem{T},r::Int64=0)

Compute the embeding of Π={bc | b ∈ k[x]/(P) , c ∈ k[y]/(Q)} ⊂ k[x,y]/(P,Q) in k[z]/(R).

# Correctness
* The algorithm is linear
"""
function embed{T}(b::Array{T,1},P::PolyElem{T},c::Array{T,1},Q::PolyElem{T},r::Int64=0)
  if r==0
    r=length(b)*length(c)
  end
  t::Array{T,1}=remTnaif(b,P,r)
  u::Array{T,1}=remTnaif(c,Q,r)
  return T[t[j]*u[j] for j in 1:r]
end

export berlekampMassey

"""
    berlekampMassey{T <: FieldElem}(a::Array{T,1},n::Int64,S=0)

Compute the minimal polynomial of a linear recurring sequence.
"""
function berlekampMassey{T <: FieldElem}(a::Array{T,1},n::Int64,S=0)
  if S==0
    k::Nemo.Field=parent(a[1])
    S::Nemo.Ring,x::PolyElem{T}=PolynomialRing(k,"x")
  else
    x=gen(S)
  end
  m::Int64=2*n-1
  R0::PolyElem{T}=S(x^(2*n))
  R1::PolyElem{T}=S(reverse(a))
  V0::PolyElem{T}=S(0)
  V1::PolyElem{T}=S(1)
  while n<=degree(R1)
    Q::PolyElem{T},R::PolyElem{T}=divrem(R0,R1)
    V::PolyElem{T}=V0-Q*V1
    V0=V1
    V1=V
    R0=R1
    R1=R
  end
  return V1*lead(V1)^(-1)
end

export computeR

"""
    computeR{T}(P::PolyElem{T},Q::PolyElem{T})

Compute the composed product R of P and Q.

# Correctness
* R is irreducible

"""
function computeR{T}(P::PolyElem{T},Q::PolyElem{T})
  m::Int64=degree(P)
  n::Int64=degree(Q)
  k::Nemo.Field=base_ring(P)

  vp::Array{T,1}=Array(T,m) # creation of the vector (1,0,...,0) (length m)
  vp[1]=k(1)
  for j in 2:m
    vp[j]=k(0)
  end

  vq::Array{T,1}=Array(T,n) # creation of the vector (1,0,...,0) (length n)
  vq[1]=k(1)
  for j in 2:n
    vq[j]=k(0)
  end

  up::Array{T,1}=monomialToDual(vp,P)
  uq::Array{T,1}=monomialToDual(vq,Q)

  t::Array{T,1}=embed(up,P,uq,Q,2*m*n)

  return berlekampMassey(t,m*n,parent(P))
end

export project

"""
    project{T}(a::Array{T,1},P::PolyElem{T},Q::PolyElem{T})

Compute the section of the embedding k[x]/(P) ⟶ k[z]/(R), where R = P ⊙ Q.

"""
function project{T}(a::Array{T,1},P::PolyElem{T},Q::PolyElem{T})
  n::Int64=degree(Q)
  m::Int64=degree(P)
  c::Array{T,1}=Array(T,n)
  k::Nemo.Field=base_ring(Q)
  c[1]=k(1)
  for j in 2:n
    c[j]=k(0)
  end
  u::Array{T,1}=remTnaif(c,Q,m*n)
  K::Nemo.Ring=parent(P)
  d::PolyElem{T}=K([a[j]*u[j] for j in 1:(m*n)])%P
  return T[coeff(d,j) for j in 0:(m-1)]
end

export phi1

"""
    phi1{T}(b::Array{T,2},P::PolyElem{T},Q::PolyElem{T})

Compute the isomorphism Φ : k[x,y]/(P,Q) ⟶ k[z]/(R).

# Arguments
* b::Array{T,2} is a 2d array reprenting (b_ji) where 0 <= i < m and 0 <= j < n, it is
the transpose of the b_ij of the paper because we prefer to extract columns than lines for
types reasons.
"""
function phi1{T}(b::Array{T,2},P::PolyElem{T},Q::PolyElem{T})
  k::Nemo.Field=parent(b[1,1])
  @assert k==base_ring(Q)
  m::Int64=degree(P)
  n::Int64=degree(Q)
  u::Array{T,1}=remTnaif(monomialToDual(T[k(1)],P),P,m*(n+1)-1)
  a::Array{T,1}=Array(T,m*n)
  for j in 1:m*n
    a[j]=k(0)
  end
  for i in 1:m
    t::Array{T,1}=remTnaif(b[:,i],Q,m*n)
    for j in 1:m*n
      a[j]=a[j]+t[j]*u[i+j-1] # /!\ indices
    end
  end
  return a
end

export inversePhi1

"""
    inversePhi1{T}(a::Array{T,1},P::PolyElem{T},Q::PolyElem{T})

Compute Φ^(-1) : k[z]/(R) ⟶ k[x,y]/(P,Q).

# Output
* b::Array{T,2} is a 2d array reprenting (b_ji) where 0 <= i < m and 0 <= j < n, it is
the transpose of the b_ij of the paper because we prefer to extract columns than lines for
types reasons.
"""
function inversePhi1{T}(a::Array{T,1},P::PolyElem{T},Q::PolyElem{T})
  k::Nemo.Field=parent(a[1,1])
  @assert k==base_ring(Q)
  K::Nemo.Ring=parent(Q)
  y::PolyElem{T}=gen(K)
  m::Int64=degree(P)
  n::Int64=degree(Q)
  b::Array{T,2}=Array(T,(n,m))
  u::Array{T,1}=remTnaif(monomialToDual(T[k(1)],P),P,m*(n+1)-1)
  for i in m:-1:1
    d::PolyElem{T}=K([a[j]*u[i+j-1] for j in 1:(m*n)])%Q
    for j in 1:n
    b[j,i]=coeff(d,j-1)
    end
  end
  return b
end

export phi2

"""
    phi2{T}(b::Array{T,2},P::PolyElem{T},Q::PolyElem{T})

Compute the isomorphism Φ : k[x,y]/(P,Q) ⟶ k[z]/(R).

# Argument
* This time b::Array{T,2} is the same as the one in the text.
"""

function phi2{T}(b::Array{T,2},P::PolyElem{T},Q::PolyElem{T})
  k::Nemo.Field=parent(b[1,1])
  @assert k==base_ring(P)
  K::Nemo.Ring=parent(P)
  z=gen(K)
  m::Int64=degree(P)
  n::Int64=degree(Q)
  N::Int64=n+m-1
  p::Int64=ceil(sqrt(N))
  q::Int64=ceil(N/p)
  y::Array{T,1}=monomialToDual(T[k(0),k(1)],Q)
  up::Array{T,1}=monomialToDual(T[k(1)],P)
  R::PolyElem{T}=computeR(P,Q)
  S::PolyElem{T}=K(dualToMonomial(embed(up,P,y,Q),R))
  U::PolyElem{T}=gcdinv(S,R)[2]

  Sprime::Array{PolyElem{T},1}=Array(PolyElem{T},q+1)

  for i in 1:(q+1)
    Sprime[i]=(S^(i-1))%R
  end

  MT::Nemo.Ring=MatrixSpace(K,q,n)
  mt::MatElem=MT()

  for i in 1:q # access matrices in column is better
    c::Array{T,1}=T[coeff(Sprime[i],h) for h in 0:(m*n-1)]
    for j in 1:n
      mt[i,j]=K(c[(j-1)*m+1:j*m]) # /!\ indices
    end
  end

  MC::Nemo.Ring=MatrixSpace(K,p,q)
  mc::MatElem=MC()

  for i in 1:p # /!\ indices
    for j in 1:q
      mc[i,j]=K(T[h+i*q+j-m-2 < 1 ? k(0) : h+i*q+j-m-2 > n ? k(0) : b[h,h+i*q+j-m-2] for h in 1:m])
    end
  end

  Mv::MatElem=mc*mt

  V::Array{PolyElem{T},1}=Array(PolyElem{T},p)
  for i in 1:p
    V[i]=sum([Mv[i,j]*z^((j-1)*m) for j in 1:n])%R
  end

  a::PolyElem{T}=K()

  for i in p:-1:1
    a=(Sprime[q+1]*a+V[i])%R
  end

  a=(a*U^(m-1))%R

  return T[coeff(a,i) for i in 0:(m*n-1)]
end

export inversePhi2

function inversePhi2{T}(a::Array{T,1},P::PolyElem{T},Q::PolyElem{T})
  k::Nemo.Field=parent(a[1])
  @assert k==base_ring(P)
  K::Nemo.Ring=parent(P)
  z=gen(K)
  m::Int64=degree(P)
  n::Int64=degree(Q)
  N::Int64=n+m-1
  p::Int64=ceil(sqrt(N))
  q::Int64=ceil(N/p)
  y::Array{T,1}=monomialToDual(T[k(0),k(1)],Q)
  up::Array{T,1}=monomialToDual(T[k(1)],P)
  R::PolyElem{T}=computeR(P,Q)
  S::PolyElem{T}=K(dualToMonomial(embed(up,P,y,Q),R))
  U::PolyElem{T}=gcdinv(S,R)[2]

  Sprime::Array{PolyElem{T},1}=Array(PolyElem{T},q+1)

  for i in 1:(q+1)
    Sprime[i]=(S^(i-1))%R
  end

  MT::Nemo.Ring=MatrixSpace(K,q,n)
  mt::MatElem=MT()

  for i in 1:q # access matrices in column is better
    c::Array{T,1}=T[coeff(Sprime[i],h) for h in 0:(m*n-1)]
    for j in 1:n
      mt[i,j]=reverse(K(c[(j-1)*m+1:j*m]),m) # reverse in order to do transposed product
    end
  end

  u::PolyElem{T}=U^(m-1)%R
  du::Int64=degree(u)

  a=T[coeff(mulT(remTnaif(a,R,2*m*n-1),u,m*n-1),j) for j in 0:(m*n-1)]

  V::Array{Array{T,1},1}=Array(Array{T,1},p)

  for i in 1:p
    V[i]=remTnaif(a,R,m*n+m-1)
    a=T[coeff(mulT(remTnaif(a,R,2*m*n-1),Sprime[q+1],m*n-1),j) for j in 0:(m*n-1)]
  end

  MV::Nemo.Ring=MatrixSpace(K,p,n)
  mv::MatElem=MV()

  for i in 1:p
    for j in 1:n
      mv[i,j]=K(V[i][(j-1)*m+1:(j-1)*m-2*m-1])
    end
  end

  mc::MatElem=mv*mt


end

end
