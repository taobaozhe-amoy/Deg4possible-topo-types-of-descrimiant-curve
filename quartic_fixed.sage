#INPUT: A homogeneous quartic in QQ[x,y,z]
#OUTPUT: The type of the real variety (from Table 1) and the 28 bitangents (with floating point coefficients)

import sys

F=QQ

T.<x,y,z>=PolynomialRing(F)

#INPUT:
f=10*x^4+15*x^3*y-17*x^2*y^2+15*x*y^3+10*y^4+15*x^3*z-833*x^2*y*z-833*x*y^2*z+15*y^3*z-17*x^2*z^2-833*x*y*z^2-17*y^2*z^2 + 15*x*z^3 + 15*y*z^3 + 10*z^4


#Uncomment the following code to produce a random quartic:
#limit=100
#Coeff=matrix([(randint(0,limit)/random_prime(limit)) for i in [0..14]])
#mon=matrix([x^4,x^3*y,x^3*z,x^2*y^2,x^2*y*z,x^2*z^2,x*y^3,x*y^2*z,x*y*z^2,x*z^3,y^4,y^3*z,y^2*z^2,y*z^3,z^4])
#f=(Coeff*transpose(mon))[0,0]

print(f)

#check whether this quartic is non-singular
Grad=ideal(f,diff(f,x),diff(f,y),diff(f,z))
if Grad.dimension()==0:
    print("Yes, it is smooth.")
else:
    sys.exit("Quartic is not smooth!")


R.<x,y,z,a,b,a0,a1,a2,a3,a4>=PolynomialRing(F)
f0=f.base_extend(R)
S.<a,b>=PolynomialRing(F)
digits=50
threshold=0.000000000001
almostzero=threshold

Line= a*x+b*y+z;
puresquare=ideal(a0*a3^2-a1^2*a4,8*a0^2*a3-4*a0*a1*a2+a1^3,8*a1*a4^2-4*a2*a3*a4+a3^3,8*a0*a1*a4-4*a0*a2*a3+a1^2*a3,8*a0*a3*a4-4*a1*a2*a4+a1*a3^2,16*a0^2*a4+2*a0*a1*a3-4*a0*a2^2+a1^2*a2,16*a0*a4^2+2*a1*a3*a4-4*a2^2*a4+a2*a3^2);
Res=f0.resultant(Line,z)
Res=Res.subs(y=1)
phi=hom(R,S,[0,0,0,a,b,Res.coefficient({x:0}),Res.coefficient({x:1}),Res.coefficient({x:2}),Res.coefficient({x:3}),Res.coefficient({x:4})])
bit1 = phi(puresquare)

I=singular.groebner(singular(bit1))
singular.lib('solve.lib')
VRing=singular.solve(I,digits)
singular.set_ring(VRing)
B1=singular("SOL")

nreal1=0
Bitangents=[]
RealBitangents=[]
for k in [1..len(B1)]:
    real=0;
    if ((B1[k][1].impart()).absValue()<threshold) and ((B1[k][2].impart()).absValue()<threshold):
        real=1
        RealBitangents=RealBitangents+[(float(B1[k][1].repart())+float(B1[k][1].impart())*i)*x+(float(B1[k][2].repart())+float(B1[k][2].impart())*i)*y+z]
    nreal1=nreal1+real
    Bitangents=Bitangents+[(float(B1[k][1].repart())+float(B1[k][1].impart())*i)*x+(float(B1[k][2].repart())+float(B1[k][2].impart())*i)*y+z]
Line=a*x+y
Res=f0.resultant(Line,y)
Res=Res.subs(z=1)
phi=hom(R,S,[0,0,0,a,0,Res.coefficient({x:0}),Res.coefficient({x:1}),Res.coefficient({x:2}),Res.coefficient({x:3}),Res.coefficient({x:4})])
bit2=phi(puresquare)+ideal(b)

if dimension(bit2)==-1:
    nreal2=0
else:
    I=singular.groebner(singular(bit2))
    singular.lib('solve.lib')
    VRing=singular.solve(I,digits)
    singular.set_ring(VRing)
    B2=singular("SOL")
    nreal2=0
    for k in [1..len(B2)]:
        real=0;
        if ((B2[k][1].impart()).absValue()<threshold) and ((B2[k][2].impart()).absValue()<threshold):
            real=1
            RealBitangents=RealBitangents+[(float(B2[k][1].repart())+float(B2[k][1].impart())*i)*x+y]
        nreal2=nreal2+real
        Bitangents=Bitangents+[(float(B2[k][1].repart())+float(B2[k][1].impart())*i)*x+y]

Res=f0.resultant(x)
Res=Res.subs(z=1)
phi=hom(R,F,[0,0,0,0,0,Res.coefficient({y:0}),Res.coefficient({y:1}),Res.coefficient({y:2}),Res.coefficient({y:3}),Res.coefficient({y:4})])
bit3=phi(puresquare)
if bit3==ideal(0):
    nreal3=1
    Bitangents=Bitangents+[x]
    RealBitangents=RealBitangents+[x]
else:
    nreal3=0

NRealBit=nreal1+nreal2+nreal3
if len(Bitangents)!=28:
    print("Something has gone wrong. We found "+str(len(Bitangents))+" bitangents")
print("The quartic has 28 bitangets, stored in 'Bitangents', and "+str(NRealBit)+" real bitangents, stored in 'RealBitangents'.")

Type=0
if NRealBit==28:
    Type="consists of 4 ovals."
if NRealBit==16:
    Type="consists of 3 ovals."
if NRealBit==8:
    Type="consists of two non-nested ovals.";

M4Bit=0
if not Type==0:
    print ("The real variety of the quartic " + Type)
    WhatDoTheyWant=input("Do you want to compute the Gram matrices of f? (y/n): ")
    if WhatDoTheyWant=="n":
        sys.exit()
    if WhatDoTheyWant=="y":
        M4Bit=1

kk=QQ
R.<a11,a12,a13,a14,a15,a16,a22,a23,a24,a25,a26,a33,a34,a35,a36,a44,a45,a46,a55,a56,a66>=PolynomialRing(kk)

A = matrix([[a11,a12,a13,a14,a15,a16],[a12,a22,a23,a24,a25,a26],[a13,a23,a33,a34,a35,a36],[a14,a24,a34,a44,a45,a46],[a15,a25,a35,a45,a55,a56],[a16,a26,a36,a46,a56,a66]])

Minors=ideal(A.minors(4))
S.<x,y,z>=PolynomialRing(R)
f=f.base_extend(S)

monomial=matrix([x^2,y^2,z^2,x*y,x*z,y*z])

Result = (monomial*A*transpose(monomial))[0,0]

GramFiber = ideal(Result.monomial_coefficient(x^4)-f.monomial_coefficient(x^4),
Result.monomial_coefficient(x^3*y)-f.monomial_coefficient(x^3*y),
Result.monomial_coefficient(x^3*z)-f.monomial_coefficient(x^3*z),
Result.monomial_coefficient(x^2*y^2)-f.monomial_coefficient(x^2*y^2),
Result.monomial_coefficient(x^2*y*z)-f.monomial_coefficient(x^2*y*z),
Result.monomial_coefficient(x^2*z^2)-f.monomial_coefficient(x^2*z^2),
Result.monomial_coefficient(x*y^3)-f.monomial_coefficient(x*y^3),
Result.monomial_coefficient(x*y^2*z)-f.monomial_coefficient(x*y^2*z),
Result.monomial_coefficient(x*y*z^2)-f.monomial_coefficient(x*y*z^2),
Result.monomial_coefficient(x*z^3)-f.monomial_coefficient(x*z^3),
Result.monomial_coefficient(y^4)-f.monomial_coefficient(y^4),
Result.monomial_coefficient(y^3*z)-f.monomial_coefficient(y^3*z),
Result.monomial_coefficient(y^2*z^2)-f.monomial_coefficient(y^2*z^2),
Result.monomial_coefficient(y*z^3)-f.monomial_coefficient(y*z^3),
Result.monomial_coefficient(z^4)-f.monomial_coefficient(z^4)
)

#This is the ideal of all rank 3 Gram matrices
VerticesIdeal=GramFiber+Minors

if VerticesIdeal.dimension()!=0:
    print("FAILURE: Gram matrix ideal has the wrong dimension!")

#Use SINGULAR to find the points
I=singular.groebner(singular(VerticesIdeal))
singular.lib('solve.lib')
VRing=singular.solve(I,digits)
singular.set_ring(VRing)
Vsing=singular("SOL")

#There should be 63 o' them

print("We have found:")
if len(Vsing)==63:
    print("63 complex Gram matrices of rank 3. Stored in 'V'.")
else:
    print("FAILURE: We did not find 63 complex Gram matrices of rank 3")

#Copy this back into Sage and into matrix-format
IndexList=[[0,0],[0,1],[0,2],[0,3],[0,4],[0,5],[1,1],[1,2],[1,3],[1,4],[1,5],[2,2],[2,3],[2,4],[2,5],[3,3],[3,4],[3,5],[4,4],[4,5],[5,5]]

V=[]
Vre=[]
Vim=[]

Zero6=copy(zero_matrix(CDF,6))

for r in [1..63]:
    k=0
    A=copy(Zero6)
    Are=copy(Zero6)
    Aim=copy(Zero6)
    for I in IndexList:
        k=k+1
        a=(Vsing[r][k]).repart()
        b=(Vsing[r][k]).impart()
        af=float(a)
        bf=float(b)
        A[I[0],I[1]]=af+bf*i
        A[I[1],I[0]]=af+bf*i
        Are[I[0],I[1]]=af
        Are[I[1],I[0]]=af
        Aim[I[0],I[1]]=bf
        Aim[I[1],I[0]]=bf
    V.append(A)
    Vre.append(Are)
    Vim.append(Aim)

Vreal=[]
Ireal=[]

for r in [0..62]:
    add=true
    for I in IndexList:
        if abs(Vim[r][I[0],I[1]])>almostzero:
            add=false
    if add:
        Vreal.append(V[r])
        Ireal.append(r)

Vpsd=[]
Ipsd=[]

for r in [0..len(Vreal)-1]:
    A=Vreal[r]
    e=A.eigenvalues()
    add=true
    for x in e:
        if x<-almostzero:
            add=false
    if add:
        Vpsd.append(A)
        Ipsd.append(r)

Vnsd=[]
Insd=[]

for r in [0..len(Vreal)-1]:
    A=Vreal[r]
    e=A.eigenvalues()
    add=true
    for x in e:
        if x>almostzero:
            add=false
    if add:
        Vnsd.append(A)
        Insd.append(r)

print(str(len(Vreal))+" real Gram matrices. Stored in 'Vreal'.")
print(str(len(Vpsd))+" psd Gram matrices. Stored in 'Vpsd'.")
print(str(len(Vnsd))+" nsd Gram matrices. Stored in 'Vnsd'.")

if len(Vreal)==7:
    Type="consists of 1 oval."
if len(Vpsd)==8:
    Type="is empty."
if len(Vnsd)==8:
    Type="is empty."
if (len(Vreal)==15) and (len(Vpsd)==0) and (len(Vnsd)==0):
    Type="consists of 2 nested ovals (a Vinnikov curve)."
if (Type==0) and (M4Bit==0):
    print("Something went wrong. No cases were satisfied. Try increasing accuracy threshold.")
else:
    if M4Bit==0:
        print ("The real variety of the quartic " + Type)
