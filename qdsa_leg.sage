reset()
var('xp,xq,mu')

xsxd = ((xp^2+xq^2+xp*xq) - (mu+1)*(xp+xq)+mu)^2 + 2*((mu+1)-(xp+xq))*((xp^3+xq^3) - (mu+1)*(xp^2+xq^2)+mu*(xp+xq)) + ((mu+1)-(xp+xq))^2*(xp-xq)^2

print xsxd.simplify_full()

