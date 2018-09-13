KLine := recformat<asq,bsq,Asq,Bsq>;

function KummerLine(asq,bsq)
	K := rec<KLine | asq := asq, bsq := bsq>;
	Asq := asq + bsq;
	Bsq := asq - bsq;
	K`Asq := Asq;
	K`Bsq := Bsq;
	return K;
end function;

function MultKL(p, k, K)
	n := Reverse(Intseq(k,2)); //print n;

	asq := K`asq; bsq := K`bsq; Asq := K`Asq; Bsq := K`Bsq;
	xsq := p[1]; zsq := p[2];

	x1 := xsq; z1 := zsq;
	x0 := Bsq*(x1+z1)^2;	//print "x0 ",(x1+z1);
	z0 := Asq*(x1-z1)^2;	//print "z0 ",(x1-z1);
	x2 := bsq*(x0+z0)^2;	//print "x2 ",x2;
	z2 := asq*(x0-z0)^2;	//print "z2 ",z2;
	i := 2;
	while IsDefined(n,i)  do
		if n[i] eq 0 then
			x0 := Bsq*(x1+z1)^2;
			z0 := Asq*(x1-z1)^2;
			x3 := bsq*(x0+z0)^2;
			z3 := asq*(x0-z0)^2;
			x0 := Bsq*(x1+z1)*(x2+z2);
			z0 := Asq*(x1-z1)*(x2-z2);
			x4 := zsq*(x0+z0)^2;
			z4 := xsq*(x0-z0)^2;
			x1 := x3; z1 := z3; x2 := x4; z2 := z4; 
		end if;
		if n[i] eq 1 then
			x0 := Bsq*(x1+z1)*(x2+z2);
			z0 := Asq*(x1-z1)*(x2-z2);
			x3 := zsq*(x0+z0)^2;
			z3 := xsq*(x0-z0)^2;
			x0 := Bsq*(x2+z2)^2;
			z0 := Asq*(x2-z2)^2;
			x4 := bsq*(x0+z0)^2;
			z4 := asq*(x0-z0)^2;
			x1 := x3; z1 := z3; x2 := x4; z2 := z4; 
		end if;
		i:=i+1;
	end while;
	return [x1,z1];
end function;

function MapToKL(P, K)
	asq := K`asq; bsq := K`bsq;
	X := P[1]/P[3];
	xsq := X*bsq;
	zsq := (X-1)*asq;
	return [xsq,zsq];
end function;

function MapToEC(p, K)
	asq := K`asq; bsq := K`bsq; 
	xsq := p[1]; zsq := p[2];
	l := asq^2/(asq^2-bsq^2);
	E := EllipticCurve([0,-l-1,0,l,0]);
	den := asq*xsq - bsq*zsq;
	if den eq 0 then return E!0; end if;
	num := asq*xsq;
	X := num/den;
	ptw := Points(E,X);
	if IsEmpty(ptw) eq false then return E!ptw[1]; end if;
end function;


 prime := 2^(251)-9; Fp := GF(prime); asq := Fp!81; bsq := Fp!20;  basept := [64, 1]; ord := 452312848583266388373324160190187140049000320168872127505022858504236695257;
// prime := 2^(255)-19; Fp := GF(prime); asq := Fp!82; bsq := Fp!77;  basept := [31, 1]; ord := 4824670384888174809315457708695329493883939011885747436657444590489242149187;
//prime := 2^(266)-3; Fp := GF(prime); asq := Fp!260; bsq := Fp!139; basept := [2, 1]; ord := 9880924948250982009478057387408034803478684522013484352184368596384732719002519;
fn := GF(ord);

K := KummerLine(asq, bsq);
l := asq^2/(asq^2-bsq^2); // print "lambda = ", l;
E := EllipticCurve([Fp| 0, -l-1, 0, l, 0]);
T2 := E![l,0,1];

//P1 := Random(E); P := P1 + T2; p := MapToKL(P, K);

// k := 3;
//k := Random(1,prime);
//print "k = ", k;

//q := MultKL(p,k,K);
//Q := k*P1 + T2;
//qprime := MapToKL(Q, K);
//q[1]/q[2]; qprime[1]/qprime[2]; // verifies KL ladder starting from a random point on the EC

//basePt := MapToEC(basept, K) + T2;
//ord*basePt + T2; // verifies order of the base point on the KL


//KeyGen
d1:=[93,97,228,133,50,83,193,84,41,86,193,35,42,170,
207,160,41,45,135,122,67,84,155,58,65,206,139,143,121,202,145];

t1:= d1[1];
for i:=2 to 31 do
	t1 := t1 + d1[i]*256^(i-1);
end for;
d1 := t1;
print "d1 := ",d1;

pk := MultKL(basept,d1,K);
pk1 := pk[1]/pk[2];
print "pk_m:=",pk1;

pk2:=[60,117,207,78,197,204,8,179,99,252,213,206,49,122,208,179,
106,188,145,4,154,6,102,111,117,111,85,79,121,37,172,3];
t1:= pk2[1];
for i:=2 to 32 do
	t1 := t1 + pk2[i]*256^(i-1);
end for;
pk2 := t1;
print "pk_c := ",pk2;
print "pk_c-pk_m:=",(pk1-pk2);

rtemp:=[123,143,49,153,102,213,232,113,106,141,234,187,252,185,124,97,3,
211,84,153,126,65,19,165,120,50,193,235,37,20,93,112,78,26,100,45,
4,12,197,68,109,74,254,134,165,190,87,45,162,30,250,178,218,68,48,
91,60,60,23,249,141,253];
t1:= rtemp[1];
for i:=2 to 62 do
	t1 := t1 + rtemp[i]*256^(i-1);
end for;
rtemp := fn!t1;
print "rtemp:=",rtemp;


r:=[195,82,71,209,122,155,42,151,247,32,39,50,225,78,169,236,241,155,
4,244,148,27,246,61,68,5,183,125,126,38,248];
t1:= r[1];
for i:=2 to 31 do
	t1 := t1 + r[i]*256^(i-1);
end for;
r := t1;
print "r := ",r;

print "r-rtemp := ",r-rtemp;

R:=[145,179,194,1,78,53,165,194,118,118,68,222,139,95,120,240,40,206,
246,206,253,98,141,156,24,164,242,241,111,7,224,4];
t1:= R[1];
for i:=2 to 32 do
	t1 := t1 + R[i]*256^(i-1);
end for;
R := t1;
print "R := ",R;

R1 := MultKL(basept,r,K);
R_m := R1[1]/R1[2];
print "R_m:=",R_m;

print "R_c-R_m := ",R-R_m;

hb:=[158,52,144,42,54,28,242,12,75,233,167,66,189,86,200,8,210,93,99,
173,67,155,107,163,1,219,19,190,146,187,218,137,75,169,120,197,110,
101,222,58,114,7,152,119,216,185,146,122,38,166,51,129,166,14,200,
14,39,130,18,22,70,82];
t1:= hb[1];
for i:=2 to 62 do
	t1 := t1 + hb[i]*256^(i-1);
end for;
hb := fn!t1;
print "hb:=",hb;

h:=[95,232,215,88,238,160,102,38,81,93,170,66,170,97,241,198,75,
173,234,81,108,43,136,82,184,155,122,12,20,67,54];
t1:= h[1];
for i:=2 to 31 do
	t1 := t1 + h[i]*256^(i-1);
end for;
h := t1;
print "h := ",h;

print "hb-h := ",hb-h;

s:=[8,198,191,178,147,37,149,211,227,50,114,157,255,254,61,4,136,
222,202,98,15,182,107,135,82,57,205,144,55,229,127];
t1:= s[1];
for i:=2 to 31 do
	t1 := t1 + s[i]*256^(i-1);
end for;
s := t1;
print "s := ",s;

s_m:= fn!(r-h*d1);
print "s_m := ",s_m;
print "s-s_m := ",s-s_m;

b51:=2^51;
sRx:=[2193493405200540,1290121454737215,214467156829684,
406732738141978,85745208440328];
sRx:=Fp!(sRx[1]+sRx[2]*b51^1+sRx[3]*b51^2+sRx[4]*b51^3+sRx[5]*b51^4);

sRz:=[1786002994812091,2205864804953583,819105131127186,
758472751897370,29882573389147];
sRz:=Fp!(sRz[1]+sRz[2]*b51^1+sRz[3]*b51^2+sRz[4]*b51^3+sRz[5]*b51^4);
sR:= sRx/sRz;
print "sR:=", sR;

sR1 := MultKL(basept,s,K);
sR_m := sR1[1]/sR1[2];
print "sR_m:=",sR_m;

print "sR-sR_m := ",sR-sR_m;

sQx:=[973599921197053,1881130351416889,1834815713981721,869388760580629,71522967936218];
sQx:=Fp!(sQx[1]+sQx[2]*b51^1+sQx[3]*b51^2+sQx[4]*b51^3+sQx[5]*b51^4);

sQz:=[877209215092906,1329101752153665,178832814841359,1396316076256463,14722978532914];
sQz:=Fp!(sQz[1]+sQz[2]*b51^1+sQz[3]*b51^2+sQz[4]*b51^3+sQz[5]*b51^4);
sQ:= sQx/sQz;
print "sQ:=", sQ;


Pt:= [Fp!pk1,1];
sQ1 := MultKL(Pt,h,K);
sQ_m := sQ1[1]/sQ1[2];
print "sQ_m:=",sQ_m;
print "sQ-sQ_m := ",sQ-sQ_m;

xpl:=[551054503923652,782863920105695,1045321024429953,637697934185303,27951142480358];

zpl:=[1717182673617704,2002032589584351,1259577515206757,1510242426086733,1896833937730];

xpl:=Fp!(xpl[1]+xpl[2]*b51^1+xpl[3]*b51^2+xpl[4]*b51^3+xpl[5]*b51^4);
zpl:=Fp!(zpl[1]+zpl[2]*b51^1+zpl[3]*b51^2+zpl[4]*b51^3+zpl[5]*b51^4);
//pl:= xpl/zpl;
//print "pl:=", pl;
print "xpl:=", xpl;
print "zpl:=", zpl;

mu:=asq^2/(asq^2-bsq^2);
xpl2:=sRz*asq;
zpl2:=xpl2-sRx*bsq;
//xpl:= xpl2/zpl2;
//ypl:= Sqrt(xpl*(xpl-1)*(xpl-mu));
//pt1:=E![xpl,ypl];
//print "pl:=", xpl2/zpl2;
//print "xpl2:=", xpl2;
//print "zpl2:=", zpl2;

xql2:=sQz*asq;
zql2:=xql2-sQx*bsq;
//xql:= xql2/zql2;
//yql:= Sqrt(xql*(xql-1)*(xql-mu));
//pt2:=E![xql,yql];

//print "pt1 := ",pt1;
//print "pt2 := ",pt2;
//print "pt1+pt2 := ",pt1+pt2;
//print "pt1-pt2 := ",pt1-pt2;

xrl2:=1*asq;
zrl2:=xrl2-R*bsq;
print "ptr := ",xrl2/zrl2;


x1:=xpl2; z1:=zpl2;
x2:=xql2; z2:=zql2;
xr:=xrl2; zr:=zrl2;
xp:=x1;zp:=z1;
xq:=x2;zq:=z2;
mu:=asq^2/(asq^2-bsq^2);


bxx:= (xp*xq - mu*zp*zq)^2;
print "bxx:= ",bxx;
bxx1:=[2153033131465093,134618518085973,283552920567959,424321949703993,90767134423769];
bxx1:=Fp!(bxx1[1]+bxx1[2]*b51^1+bxx1[3]*b51^2+bxx1[4]*b51^3+bxx1[5]*b51^4);
print "bxx1:= ",bxx1;

bxz := (xp*zq+xq*zp)*(xp*xq+mu*zp*zq)-2*(mu+1)*(xp*xq)*zp*zq;
print "bxz:= ",bxz;
bxz1:=[2099432947353384,2096864751155795,1246802994607786,1870683091594301,131291624550327];
bxz1:=Fp!(bxz1[1]+bxz1[2]*b51^1+bxz1[3]*b51^2+bxz1[4]*b51^3+bxz1[5]*b51^4);
print "bxz1:= ",bxz1;


bzz := (xp*zq-zp*xq)^2;
print "bzz:= ",bzz;
bzz1:=[1093173696694666,813838379638443,1819301496764920,1691003426661428,72420303276951];
bzz1:=Fp!(bzz1[1]+bzz1[2]*b51^1+bzz1[3]*b51^2+bzz1[4]*b51^3+bzz1[5]*b51^4);
print "bzz1:= ",bzz1;


f := bzz*xr^2-2*bxz*xr*zr+bxx*zr^2;
ZZ:=IntegerRing();
print "f:= ",f;//Intseq(ZZ!f,2^51);

//print Intseq(ZZ!mu,2^51):Hex;
print Intseq(prime,2^51);


f1:=[804485027731767,  1939189742132221,  644969960674411,  858583021524075,  106493690151423];
f1:=(f1[1]+f1[2]*b51^1+f1[3]*b51^2+f1[4]*b51^3+f1[5]*b51^4);
print "f1:= ",Fp!f1;
print "f1-p:= ",f1-prime;

print (2^52-1);








