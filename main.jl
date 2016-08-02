using Nemo

function monomialToDual(a,P)
  """
  Convert monomial (canonical) coordinates "a" to dual coordinates
  with respect to the trace engendered by P.

    INPUT:
  - a : list of m elements in k (perfect field) ;
  - P : polynom of degree m over k, P must be monic and squarefree.

    OUTPUT:
  - list of m elements in k (which represent the coordinates in the dual basis).

    REMARK:
  k[x]/(P) is not necessarly a field since P is not necessarly irreducible.
  """
  k=parent(a[1])
  m=degree(P)
  R=parent(P)
  t=gen(R)
  # R,t=PolynomialRing(k,"t") # creation of polynom space and power series space
  S,x=PowerSeriesRing(k,m,"x")
  Q=reverse(P,m+1)
  Q=S([k(coeff(P,i)) for i in 0:(m-1)],m,m) # the k(coeff(...)) is here to be sure that we work with the correct type
  T=inv(Q) # here Q belongs to PowerSeriesRing and is inversible
  T=R([k(coeff(T,i)) for i in 0:(m-1)]) # here T goes from PowerSeriesRing to PolynomialRing
  A=R(a)
  b=(reverse((derivative(P)*A)%P,m)*T)%(t^m)
  return [k(coeff(b,i)) for i in 0:(m-1)]
end

function dualToMonomial(b,P)
  """
  Convert dual coordinates "b" to monomial coordinates
  with respect to the trace engendered by P.

    INPUT:
  - b : list of m elements in k (perfect field) ;
  - P : polynom of degree m over k, P must be monic and squarefree.

    OUTPUT:
  - list of m elements in k (which represent the coordinates in the monomial basis).

    REMARK:
  k[x]/(P) is not necessarly a field since P is not necessarly irreducible.
  """
  k=parent(b[1])
  m=degree(P)
  R=parent(P)
  t=gen(R)
#  R,t=PolynomialRing(k,"t")
  S=gcdinv(derivative(P),P)[2]
  c=(reverse(P,m+1)*R(b))%(t^m)
  c=reverse(c,m)
  d=(c*S)%P
  return [k(coeff(d,i)) for i in 0:(m-1)]
end

function mulT(c,P,n)
  """
  The tranposition of the algorithm of multiplication by P.

    INPUT:
  - c : a list of m+n elements in k (perfect field)
  - P : a polynom over k
  - n : integer, needed ?
  """
  m=degree(P)
  k=parent(coeff(P,0))
  R,t=PolynomialRing(k,"t")
  p=fq_nmod[k(coeff(P,j)) for j in 0:m]
  b=Array(fq_nmod,(1,n+1)) # un tableau avec des #undef dedans, la fonction zeros ne fonctionne pas avec fq_nmod
  for i in 1:(n+1)
    b[i]=k(0) # comme la fonction zeros ne fonctionne, on met à 0 à la main
  end
  for i in range(m+n,-1,m+n+1)
    for j in range(min(m,i),-1,min(m,i)-max(0,i-n)+1)
      println((i,j))
      b[i-j+1]=b[i-j+1]+p[j+1]*c[i+1]
    end
  end
  return R(fq_nmod[k(b[i]) for i in 1:(n+1)])
end

function remT(r,P,n)
  """
    Ne fonctionne pas !! À débugger.
  """
  m=degree(P)
  k=parent(coeff(P,0))
  T,t=PolynomialRing(k,"t") # creation of polynom space and power series space
  S,x=PowerSeriesRing(k,m,"x")
  R=T(fq_nmod[r[i] for i in 1:m])# r est la liste des coeff et R le polynôme
  α=reverse(P,m+1)
  α=inv(S(fq_nmod[coeff(α,i) for i in 0:m],m,m))
  α=T(fq_nmod[coeff(α,i) for i in 0:(n+m)])
  b=copy(r)
  while length(b)<(n+1)
    push!(b,k(0))
  end
  return R-(t^m)*((α*mulT(b,P,n-m))%(t^(n-m+1)))
end

function remTnaif(r,P,n)
  """
  Linear extension algorithm. Takes the m first elements of a Linear
  recurring sequence and computes the n first.

    INPUT:
  - r : m first elements of a linear recurring sequence ;
  - P : the polynom of the sequence ;
  - n : integer, the number of elements desired.

    OUPUT:
  - a list of the n first elements of the linear recurring sequence.
  """
  m=degree(P)
  p=fq_nmod[coeff(P,j) for j in 0:m]
  b=copy(r)
  while length(b)<n
      s=sum([-1*p[j]*b[end-m+j] for j in 1:m])
      push!(b,p[end]^(-1)*s)
  end
  return b
end

function embed(b,P,c,Q,r=0)
  """
  Compute the embeding of Π={bc | b ∈ k[x]/(P) , c ∈ k[y]/(Q)} ⊂ k[x,y]/(P,Q)}.
  """
  if r==0
    r=length(b)*length(c)
  end
  t=remTnaif(b,P,r)
  u=remTnaif(c,Q,r)
  return fq_nmod[t[j]*u[j] for j in 1:r]
end

function berlekampMassey{T}(a::Array{T,1},n::Int64,S=0)
  if S==0
    k=parent(a[1])
    S,x=PolynomialRing(k,"x")
  else
    x=gen(S)
  end
  polyT=typeof(x)
  m::Int64=2*n-1
  R0::polyT=S(x^(2*n))
  R1::polyT=S(reverse(a))
  V0::polyT=S(0)
  V1::polyT=S(1)
  while n<=degree(R1)
    Q::polyT,R::polyT=divrem(R0,R1)
    V::polyT=V0-Q*V1
    V0=V1
    V1=V
    R0=R1
    R1=R
  end
  return V1*lead(V1)^(-1)
end

function computeR{polyT}(P::polyT,Q::polyT)
  m::Int64=degree(P)
  n::Int64=degree(Q)
  T=typeof(coeff(P,0))
  k=parent(coeff(P,0))

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

  up=monomialToDual(vp,P)
  uq=monomialToDual(vq,Q)

  t=embed(up,P,uq,Q,2*m*n)

  return berlekampMassey(t,m*n,parent(P))
end

function project(a::Array{fq_nmod,1},P::fq_nmod_poly,Q::fq_nmod_poly)
  n::Int64=degree(Q)
  m::Int64=degree(P)
  c::Array{fq_nmod,1}=Array(fq_nmod,n)
  k=parent(coeff(Q,0))
  c[1]=k(1)
  for j in 2:n
    c[j]=k(0)
  end
  u::Array{fq_nmod,1}=remTnaif(c,Q,m*n)
  println(u)
  K=parent(P)
  d=K([a[j]*u[j] for j in 1:(m*n)])%P
  return fq_nmod[coeff(d,j) for j in 0:(m-1)]
end

### ESPACE DE TESTS ###

k,u=FiniteField(5,1,"u")
T,t=PolynomialRing(k,"t")

P=t^3+t+1
Q=t^2+t+1

parent(P)

R=computeR(P,Q)

R+P

uq=monomialToDual(fq_nmod[k(1),k(0)],Q)
up=monomialToDual(fq_nmod[k(1),k(0),k(0)],P)
a=embed(up,P,uq,Q)
aa=dualToMonomial(a,R)


project(aa,P,Q)


t=Array(fq_nmod,3)
typeof(t)
t[2]=k(2)



t
