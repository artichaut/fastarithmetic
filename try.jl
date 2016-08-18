using Nemo, FastArithmetic

k,u=FiniteField(5,1,"u")
K,t=PolynomialRing(k,"t")
P=t^2+t+1
Q=t^3+t+1

b=fq_nmod[k(1) k(2);k(4) k(2); k(1) k(3)]
a=phi1(b,P,Q)

b

a2=phi2(transpose(b),P,Q)

inversePhi2(a2,P,Q)

c=fq_nmod[k(1) k(0) k(0); k(0) k(1) k(0)]

phi2(c,P,Q)

M=MatrixSpace(K,2,2)

f=M()

f[1,2]=t
f
1+1

A=Array(Array{fq_nmod,1},5)

A[1]=[k(0)]




A


for i in p:-1:1
  for j in q:-1:1
    for k in n:-1:1
      A[i,k]=A[i,k]+mulT(C[i,j],B[k,j],le bon argument)
    end
  end
end
