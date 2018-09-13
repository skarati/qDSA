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


 //prime := 2^(251)-9; Fp := GF(prime); asq := Fp!81; bsq := Fp!20;  basept := [64, 1]; ord := 452312848583266388373324160190187140049000320168872127505022858504236695257;
// prime := 2^(255)-19; Fp := GF(prime); asq := Fp!82; bsq := Fp!77;  basept := [31, 1]; ord := 4824670384888174809315457708695329493883939011885747436657444590489242149187;
prime := 2^(266)-3; Fp := GF(prime); asq := Fp!260; bsq := Fp!139; basept := [2, 1]; ord := 9880924948250982009478057387408034803478684522013484352184368596384732719002519;
fn := GF(ord);

K := KummerLine(asq, bsq);
l := asq^2/(asq^2-bsq^2); // print "lambda = ", l;
E := EllipticCurve([Fp| 0, -l-1, 0, l, 0]);
T2 := E![l,0,1];



//KeyGen
d1:=[93,97,228,133,50,83,193,84,41,86,193,35,42,170,207,160,41,45,
135,122,67,84,155,58,65,206,139,143,121,202,17,2];

t1:= d1[1];
for i:=2 to 32 do
	t1 := t1 + d1[i]*256^(i-1);
end for;
d1 := t1;
print "d1 := ",d1;

pk := MultKL(basept,d1,K);
pk1 := pk[1]/pk[2];
print "pk_m:=",pk1;

pk2:=[75,160,4,149,249,183,208,154,231,248,109,1,81,72,92,179,97,168,
0,142,182,224,197,247,101,19,124,143,234,85,7,83];
t1:= pk2[1];
for i:=2 to 32 do
	t1 := t1 + pk2[i]*256^(i-1);
end for;
pk2 := t1;
print "pk_c := ",pk2;
print "pk_c-pk_m:=",(pk1-pk2);

rtemp:=[123,143,49,153,102,213,232,113,106,141,234,187,252,185,
124,97,3,211,84,153,126,65,19,165,120,50,193,235,37,20,93,112,
78,26,100,45,4,12,197,68,109,74,254,134,165,190,87,45,162,30,
250,178,218,68,48,91,60,60,23,249,141,253,44,116];
t1:= rtemp[1];
for i:=2 to 64 do
	t1 := t1 + rtemp[i]*256^(i-1);
end for;
rtemp := fn!t1;
print "rtemp:=",rtemp;


r:=[89,161,38,165,109,64,110,66,198,37,0,78,149,201,132,34,175,190,
157,112,198,192,124,137,217,64,107,102,186,197,53,0];
t1:= r[1];
for i:=2 to 31 do
	t1 := t1 + r[i]*256^(i-1);
end for;
r := t1;
print "r := ",r;

print "r-rtemp := ",r-rtemp;

R:=[104,15,118,69,252,0,233,78,55,32,205,245,226,255,53,249,15,243,180,
203,220,109,128,108,102,189,13,132,197,65,184,38];
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

hb:=[98,192,123,93,135,142,203,197,62,37,165,194,226,120,246,53,86,
155,17,242,195,78,145,169,23,79,185,181,165,124,85,143,167,118,47,
30,181,166,226,128,104,0,181,43,139,59,49,233,37,200,42,74,75,159,
63,169,72,251,1,104,216,107,63,213];
t1:= hb[1];
for i:=2 to 64 do
	t1 := t1 + hb[i]*256^(i-1);
end for;
hb := fn!t1;
print "hb:=",hb;

h:=[9,255,50,33,33,3,137,3,116,215,227,174,75,151,58,10,68,225,
148,249,221,102,193,155,60,65,19,182,23,19,168,7];
t1:= h[1];
for i:=2 to 32 do
	t1 := t1 + h[i]*256^(i-1);
end for;
h := t1;
print "h := ",h;

print "hb-h := ",hb-h;

s:=[98,135,172,234,88,134,6,175,178,21,75,150,159,129,167,213,3,
180,222,123,32,71,126,21,67,66,135,243,148,247,15,8];
t1:= s[1];
for i:=2 to 32 do
	t1 := t1 + s[i]*256^(i-1);
end for;
s := t1;
print "s := ",s;

s_m:= fn!(r-h*d1);
print "s_m := ",s_m;
print "s-s_m := ",s-s_m;

b51:=2^54;
sRx:=[1801260713469259,65898346828219,234008582606192,898680999695052,1279602879159694];
sRx:=Fp!(sRx[1]+sRx[2]*b51^1+sRx[3]*b51^2+sRx[4]*b51^3+sRx[5]*b51^4);

sRz:=[6099672911370,1991491930751827,1620403312817891,2084496682000930,925935327156023];
sRz:=Fp!(sRz[1]+sRz[2]*b51^1+sRz[3]*b51^2+sRz[4]*b51^3+sRz[5]*b51^4);
sR:= sRx/sRz;
print "sR:=", sR;

sR1 := MultKL(basept,s,K);
sR_m := sR1[1]/sR1[2];
print "sR_m:=",sR_m;

print "sR-sR_m := ",sR-sR_m;

sQx:=[1721259115114249,758227743822392,948885451747450,27503972961140,884997816329421];
sQx:=Fp!(sQx[1]+sQx[2]*b51^1+sQx[3]*b51^2+sQx[4]*b51^3+sQx[5]*b51^4);

sQz:=[1455977701676055,1400591270033755,948388409081310,787656701462342,983865673753158];
sQz:=Fp!(sQz[1]+sQz[2]*b51^1+sQz[3]*b51^2+sQz[4]*b51^3+sQz[5]*b51^4);
sQ:= sQx/sQz;
print "sQ:=", sQ;


Pt:= [Fp!pk1,1];
sQ1 := MultKL(Pt,h,K);
sQ_m := sQ1[1]/sQ1[2];
print "sQ_m:=",sQ_m;
print "sQ-sQ_m := ",sQ-sQ_m;

xpl:=[500173178732967,1172751736311958,16882643637502,2043741897682719,1617302975180777];

zpl:=[1414686690084564,602178657909529,12620292442700,399299331721147,2167073082035221];
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
print "xpl2:=", xpl2;
print "zpl2:=", zpl2;

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
bxx1:=[1904847788134704,710625349192459,915907170910523,2187341000132357,1362838143301883];
bxx1:=Fp!(bxx1[1]+bxx1[2]*b51^1+bxx1[3]*b51^2+bxx1[4]*b51^3+bxx1[5]*b51^4);
print "bxx1:= ",bxx1;

bxz := (xp*zq+xq*zp)*(xp*xq+mu*zp*zq)-2*(mu+1)*(xp*xq)*zp*zq;
print "bxz:= ",bxz;
bxz1:=[674721113434586,1154175960894568,1810380620238074,18368970768741,2065435743030494];
bxz1:=Fp!(bxz1[1]+bxz1[2]*b51^1+bxz1[3]*b51^2+bxz1[4]*b51^3+bxz1[5]*b51^4);
print "bxz1:= ",bxz1;


bzz := (xp*zq-zp*xq)^2;
print "bzz:= ",bzz;
bzz1:=[1223705165743284,1770724106717456,120708341079672,401436238952708,504745956138034];
bzz1:=Fp!(bzz1[1]+bzz1[2]*b51^1+bzz1[3]*b51^2+bzz1[4]*b51^3+bzz1[5]*b51^4);
print "bzz1:= ",bzz1;


f := bzz*xr^2-2*bxz*xr*zr+bxx*zr^2;
ZZ:=IntegerRing();
print "f:= ",f;//Intseq(ZZ!f,2^51);

//print Intseq(ZZ!mu,2^51):Hex;
print Intseq(prime,2^51);


f1:=[18014398509481981,18014398509481983,18014398509481983,18014398509481983,1125899906842623];
f1:=(f1[1]+f1[2]*b51^1+f1[3]*b51^2+f1[4]*b51^3+f1[5]*b51^4);
print "f1:= ",Fp!f1;
print "f1-p:= ",f1-prime;

//print (2^52-1);







