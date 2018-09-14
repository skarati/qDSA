reset()

var('U1, U2, X1, Z1, X2, Z2, mu, X, Z, sumofroots, productofroots, US')

def f(U,mu):
    return U*(U-1)*(U-mu)

###############################
# Case when U1 not equal to U2#
###############################

    
def U1neqU2():
    sumofroots = 2*( mu+1 - U1 - U2 + ( f(U1, mu) + f(U2, mu) ) / (U2 - U1)^2 )
    tmp = sumofroots.simplify_full()
    tmpA = mu+1 - (U1+U2)
    V1sq = f(U1,mu)
    V2sq = f(U2,mu)
    tmpB = 2*( V2sq + V1sq )/(U2-U1)^2
    tmpC = ( V2sq - V1sq )^2 / (U2 - U1)^4
    productofroots = tmpA^2 + tmpA * tmpB + tmpC
    tmp = productofroots.simplify_full()
    tmp1 = sumofroots.substitute(U1 = X1/Z1, U2 = X2/Z2)
    tmp2 = productofroots.substitute(U1 = X1/Z1, U2 = X2/Z2)
    return [tmp1, tmp2]


###########################
# Case when U1 equal to U2#
###########################

def U1eqU2():
    US = mu+1 - 2*U1 + ( 3*U1^2 - 2*(mu+1)*U1 + mu  )^2 / ( 4*f(U1,mu) )
    tmp = US.simplify_full()
    tmp = tmp.substitute(U1 = X1/Z1, U2 = X2/Z2)
    return tmp

print "Case when U1 not equal to U2"
tmp = U1neqU2()
tmp1 = tmp[0]
tmp2 = tmp[1]
print tmp1.factor()
print "***************************"
print tmp2.factor()
tmp = tmp2

print "***************************"
print "Case when U1 equal to U2"
tmp = U1eqU2()
print "***************************"
print tmp.factor()

