import numpy as np
import math as math
import originpro as op
def RMS(tv,numv,yk):
    nn=np.size(yk)
    st=[]
    #dla każdego przedziału
    for nv in range(numv):
        #wyznaczamy dane
        xr= yk[nv*tv:((nv+1)*tv)]
        x=np.arange(nv*tv+1,(nv+1)*tv+1)
        #metoda najmniejszych kwadratów (linia zbiegająca do przedziału)
        pv=np.polyfit(x,xr,1)
        #wzór funkcji liniowej
        fv=np.poly1d(pv)
        #wyliczenie wartości funkci dla x
        ord=fv(x)
        st=np.append(st, (ord - xr)**2)
    wartosc=np.cumsum(st)
    rms=(wartosc[np.size(wartosc)-1]*1/nn)**  .5
    return rms


def DFA(xs):
    n=len(xs)
    m=np.mean(xs)     #średnia
    #suma bieżąca
    xsm=xs-m
    sums=xsm[0]
    yk=[]
    for i in range(n):
        yk.append(sums)
        sums+=xsm[i]
    sugp=[]
    for s in range(11,int(n/4)):
        if int(4*(2**(1/8))**s+0.5)<=int(n/4):
            sugp=np.append(sugp,int(4*(2**(1/8))**s+0.5))
        else:
            break

    #wartości Fn
    Fn=[]
    for tv in sugp:
        Fn=np.append(Fn,RMS(int(tv),int(n/tv),yk))
    
    prlog=np.polyfit(np.log10(sugp),np.log10(Fn),1)
    #wzór funkcji liniowej
    flog=np.poly1d(prlog)
    st_flog=str(object=flog)
    logx=np.log10(sugp)
    logy=flog(np.log10(sugp))
    logF=np.log10(Fn)
    hurst=[]
    hurst.append(round(flog[1]/2,3))
    wks = op.find_sheet('w')
    wks.del_col(1,6)
    a=wks._check_add_cols(4,2)  
    #uzupełnianie kolumn danymi
    h=wks.from_list(1,'Hurst','Hurst exponent')
    hu=wks.from_list(2,hurst)
    l=wks.from_list(3,logx,'Data scale (log)') 
    ly=wks.from_list(4,logy,'Fit')
    f=wks.from_list(5,logF,'RMS (log)') 
    gl = op.new_graph('Linear fit')[0]
    dp = gl.add_plot(wks,'F','D',type='s')
    g2 = op.new_graph('Data')[0]
    d2 = g2.add_plot(wks,'A','#') 
    d=gl.add_plot(wks,'E','D',type='l')
    d.color='#FA3C3C'
    gl.rescale()
    g2.rescale()
    return logx
    
    
def avg(paramvector):
    d=0.0
    b=0
    while(b<len(paramvector)):
        double_=paramvector[b]
        d+=double_
        b+=1
    return d/len(paramvector)

def RS(paramvector):
    d1=avg(paramvector)
    double_=paramvector[0]
    d2=float(double_)
    d3=d2
    d4=0.0
    b=0
    while (b<len(paramvector)):
           double_=paramvector[b]
           d4+=float(double_)-d1
           if (d4> d2):
               d2=d4
           elif (d4<d3):
               d3=d4
           b+=1
    return d2-d3

def var (paramvector):
    d1=avg(paramvector)
    d2=0.0
    d3=0.0
    b=0
    while (b<len(paramvector)):
        double_=paramvector[b]
        d2+=math.pow(float(double_-d1),2.0)
        b+=1
    if (len(paramvector)>1):
        d3=d2/(len(paramvector)-1)
    return d3


def RS_Plot (paramvector):
    i=int(math.floor(len(paramvector)/2))
    j=4
    vector=[]
    b=0
    while (j<=i):
        vector.append(j)
        j+=1
    vector1=[]
    
    b=0
    while (  b< len(vector)):
        k=int(vector[b])
        m=len(paramvector)/k
        b1=0
        n=0
        vector2=[]
        while(b1<math.floor(m)):
            n=b1*k
            vector3=[]
            while ((n<(b1+1)*k) and n<=len(paramvector)):
                vector3.append(paramvector[n])
                n+=1
            d1=RS(vector3)
            d2=var(vector3)
            d3=0.0
            if (d2>0.0):
                d3=d1/math.sqrt(d2)
            else:
                d3=((len(vector3)-1)/math.sqrt(len(vector3)))
            vector2.append(d3)
            b1+=1
        vector1.append(avg(vector2))
        b+=1
    b=0
    logY=[]
    logX=[]
    logF=[]
    while (b<len(vector1)):
        d1=float(vector1[b])
        k=(vector[b])
        d2=math.log(d1)
        logY.append(d2)
        d3=math.log(k)
        logX.append(d3)
        b+=1
    prlog=np.polyfit(logX,logY,1)
    #wzór funkcji liniowej
    flog=np.poly1d(prlog)
    logF=flog(logX)
    d1=0.0
    d2=0.0
    d3=0.0
    d4=0.0
    b=0
    while (b<len(logX)):
        double_1=logX[b]
        double_2=logY[b]
        d9=double_1
        d10=double_2
        d1+=d9
        d2+=d10
        d3+=d9*d9
        d4+=d9*d10
        b+=1
    d6=(d4-d1*d2/len(logX))/(d3-d1*d1/len(logX))
    h=[]
    h.append(d6)
    
    wks = op.find_sheet('w')
    wks.del_col(1,6)
    a=wks._check_add_cols(4,2)
    hh=wks.from_list(1,'Hurst','Hurst exponent')
    hu=wks.from_list(2,h)
    l=wks.from_list(3,logX,'Log(x)') 
    ly=wks.from_list(4,logY,'Log(y)')
    l2=wks.from_list(5,logF,'Fit')
    gl = op.new_graph('Linear fit')[0]
    dp = gl.add_plot(wks,'E','D',type='s')
    dd = gl.add_plot(wks,'F','D',type='l')
    dd.color='#FA3C3C'
    dp.symbol_size=3
    g2 = op.new_graph('Data')[0]
    d2 = g2.add_plot(wks,'A','#')
    gl.rescale()
    g2.rescale()
    return h