import math,sys
fn=sys.argv[1] if len(sys.argv)>1 else "vortex_horn_torus.stl"
tris=[];cur=[]
for line in open(fn):
    s=line.split()
    if s and s[0]=="vertex":
        cur.append((float(s[1]),float(s[2]),float(s[3])))
        if len(cur)==3: tris.append(tuple(cur));cur=[]
def sub(a,b):return(a[0]-b[0],a[1]-b[1],a[2]-b[2])
def cr(a,b):return(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dt(a,b):return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def nm(a):
    m=math.sqrt(dt(a,a));return (a[0]/m,a[1]/m,a[2]/m) if m else (0,0,0)
def ad(v):return math.hypot(v[0],v[2])
data=[]
for t in tris:
    a,b,c=t;n=nm(cr(sub(b,a),sub(c,a)))
    cen=((a[0]+b[0]+c[0])/3,(a[1]+b[1]+c[1])/3,(a[2]+b[2]+c[2])/3)
    data.append((a,b,c,n,cen,min(ad(v) for v in t)))
def hit(o,d,skip):
    best=1e18
    for k,(a,b,c,n,cen,md) in enumerate(data):
        if k==skip:continue
        e1=sub(b,a);e2=sub(c,a);h=cr(d,e2);det=dt(e1,h)
        if -1e-9<det<1e-9:continue
        inv=1/det;s=sub(o,a);u=dt(s,h)*inv
        if u<-1e-7 or u>1+1e-7:continue
        q=cr(s,e1);v=dt(d,q)*inv
        if v<-1e-7 or u+v>1+1e-7:continue
        tt=dt(e2,q)*inv
        if 1e-4<tt<best:best=tt
    return best
ths=[];minrr=1e9
for k,(a,b,c,n,cen,md) in enumerate(data):
    minrr=min(minrr,md)
    if md>16:continue
    o=(cen[0]-1e-4*n[0],cen[1]-1e-4*n[1],cen[2]-1e-4*n[2])
    t=hit(o,(-n[0],-n[1],-n[2]),k)
    if t<1e17: ths.append(t)
ths.sort()
print("min radius from axis (throat hole):",round(minrr,3),"mm")
print("near-throat probes",len(ths),"| min wall",round(ths[0],3),
      "| <0.8mm:",sum(1 for t in ths if t<0.8),"| <0.4mm:",sum(1 for t in ths if t<0.4))
