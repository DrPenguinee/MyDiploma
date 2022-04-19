import numpy as np
from scipy import interpolate
from scipy import integrate

import csv


# COMMON blocks:

# common/const
em = 510999.06
alpha = 1./137.0359895

# common/source
thick = .0

# common/resfun
shift =  .0
shift2 = .0

# common/lattice
epoint = np.empty(301)
npoint = 0

# common/energy_ml
pe = .0
pe0 = .0

# common/levdat
elev = np.zeros(193)
prolev = np.zeros(193)
ist = 0

# common/energy
pe = .0
pe0 = .0
psnm = .0
ntype = .0

# common/trapco
e1 = .0
e2 = .0

# common/temp
tspec = .0
scprob = .0

# common/scatt
e, s1, s2, s3, s4, s5 = np.zeros(221), np.zeros(221), np.zeros(221), np.zeros(221), np.zeros(221), np.zeros(221)

# common/theta_max
theta_max = 0.

# common/DESPL
sxt, E0 = 0., 0.


#------------------------------------------



# SAVE blocks:

# chi2-------------
ifl_chi2 = 0
Edata = np.zeros(301)
error = np.zeros(301)

# specmlint------------------
ifl_specmlint = 0
de, epoint = np.zeros(301), np.zeros(301)
slev1, cm1 = np.zeros(301), np.zeros(301)

# specint---------------------
ifl_specint = 0
thickold = .0
e0old = .0
shiftold = 100.
shift2old = 100.
elev1, prolev1 = np.zeros(193), np.zeros(193)
de, enode, epoint = np.zeros(151), np.zeros(151), np.zeros(301)
slev1, cm1 = np.zeros(151), np.zeros(151)
vsnm = np.zeros(110)
vneut1 = np.zeros(151*110).reshape((151, 110))
# cm201 = np.zeros(151*110).reshape((151, 110))
# cm221 = np.zeros(151*110).reshape((151, 110))
# cm021 = np.zeros(151*110).reshape((151, 110))

# scatprob------------------------
en, s1, g1 = np.zeros(221), np.zeros(221), np.zeros(221)
ifl_scatprob = 0

# ProbCalc_int----------------------------------
theta, P0, P1, P2, P3, P4, P5 = np.zeros(40), np.zeros(40), np.zeros(40), np.zeros(40), np.zeros(40), np.zeros(40), np.zeros(40)

#transmission--------------------------------------
thickold_t, shiftold_t = 0., 100.
spline_t = lambda x: 1




#------------------------------------------------



def chi2(xval):
    # common------------------
    global thick # source
    global shift, shift2 # resfun
    global epoint, npoint # lattice
    # ------------------------

    # save--------------------
    global ifl_chi2
    global Edata, error
    # ------------------------

    iflag = 5


    epoint, Edata, error, espnm = np.zeros(301), np.zeros(301), np.zeros(301), np.zeros(301)
    espec, backval, endef, eml = np.zeros(301), np.zeros(301), np.zeros(301), np.zeros(301)
    trspec, hnu, espnm2, fakesp = np.zeros(301), np.zeros(301), np.zeros(301), np.zeros(301)
    rancor, tspec = np.zeros(301), np.zeros(301)

    e0 = xval[0]+18570.0                # e0 - spectrum endpoint
    s = xval[1]                         # s - spectrum square
    snm = xval[2]                       # snm - squared neutrino mass
    back = xval[3]                      # back - background value
    backpar = xval[4]                   # backpar - background parameter
    step = xval[5]                      # step - value of endpoint effect
    eend = xval[6] + 18000.0            # eend - endpoint of step effect
    shift = xval[7]                     # shift - probability shift of excitation-ionization losess
    shift2 = xval[8]                    # shift2 - mainstate probability shift
    thick = xval[9]                     # thick - source thickness factor
    emin = xval[10]                     # emin - lower limit of analysis interval
    emax = xval[11]                     # emax - upper limit of analysis interval
    dtime = xval[12]                    # dtime - counting system deadtime
    dtime0 = 6.e-6
    ntype = int(xval[13] + 0.1)         # ntype = 1 - calculation by splines
                                        # ntype = 2 - direct calculation
    snm2 = xval[14]**2                  # snm2 - squared heavy neutrino mass
    hnupr = xval[15]                    # hnupr - heavy neutrino probability
    unvis = xval[16]                    # unvis - unvisible correction shift
    cut = xval[17]                      # cut - point "glitches" cut
    pos_ml = e0 - xval[18]              # pos_ml - missing level energy
    prob_ml = xval[19]                  # prob_ml - missing level probability
    acc = .0002

    if ifl_chi2 == 0:
        ifl_chi2 += 1
        
        Num_MC = 0                           # MC simulations counter

        # data file read --------------------------------------------------------
        npoint = 27

        with open('knm1.dat', 'r') as f:
            reader = csv.reader(f, delimiter=' ')
            i = 0
            for row in reader:
                epoint[i] = 1000 * float(row[0])
                Edata[i] = float(row[1])
                error[i] = float(row[2])
                print("%3d %10.1f %11.6f %13.8f" % (i, epoint[i], Edata[i], error[i]))
                i+=1
                

        
        # Arrays calculation
        for i in range(1, 101):
            entr = float(i)
            tran = transmission(entr) # need transmission -> rfcal -> lospat
            # fstail(e0 - entr - 165., e0, tail)
            vesp, vespnm = specint(e0 - entr, e0, snm, thick)
            print(entr, tran)
        # Arrays calculated

        # Faked spectrum formation--------------------------------------
        w = 96.4
        with open('run_kat_MC2.dat', 'w') as f:
            print("Simulation with Elow=%f" % (emin), file=f)
            print("%s %s %s %s" % ("HV", "Freq", "Err", "Time"), file=f)
            for i in range(npoint):
                e = epoint[i]
                espec[i], espnm[i] = specint(e, e0, snm, thick) # ok
                hnu[i], espnm2[i] = specint(e, e0, snm2, thick) # ok
                eml[i] = specmlint(e, pos_ml, e0) # ok
                # trapbackground(e, trspec[i], e0)
                endef[i] = endeffect(e, eend, step)
                backval[i] = background(e, e0, back, backpar)
                # backval[i] = backval[i]
                fakesp[i] = w * (espec[i] + hnu[i]*hnupr + eml[i]*prob_ml) + backval[i] + endef[i]
                Edata[i] = fakesp[i]
                print(e, Edata[i], error[i], file=f)
        # Faked spectrum formed--------------------------------------------

    if iflag == 5:
        Num_MC += 1
        for i in range(npoint):
            Edata[i] = fakesp[i] + error[i]*np.random.normal()
            rancor[i] = (Edata[i] - fakesp[i])/error[i]
        for i in range(npoint):
            print(rancor[i])
        print(Num_MC)

    # string 188
    if iflag == 3:
        pass

    fval = 0.
    sespec = 0.
    sml = 0.
    strap = 0.
    sbackval = 0.
    sendef = 0.
    for i in range(npoint):
        e = epoint[i]
        if iflag != 3:
            if e < emin or e > emax:
                continue
        if ntype == 1:
            espec[i], espnm[i] = specint(e, e0, snm, thick)
            hnu[i], espnm[i] = specint(e, e0, snm, thick)
        if ntype == 2:
            expspectrum(e, espec[i], e0, snm, acc) # need transmission
            expspectrum(e, hnu[i], e0, snm2, acc)
        eml[i] = specmlint(e, pos_ml, e0)
        backval[i] = background(e, e0, back, backpar)
        endef[i] = endeffect(e, eend, step) # need transmission
        # if
        if not(e < emin - .01 or e > emax):
            sespec += espec[i] + hnu[i]*hnupr
            sml += eml[i]*prob_ml
            sbackval += backval[i]
            sendef += endef[i]
    w = (s - sbackval)/sespec
    for i in range(npoint):
        e = epoint[i]
        theor = w * (espec[i] + hnu[i]*hnupr + trspec[i]*backpar + eml[i]*prob_ml) \
            + backval[i] + endef[i]
        exper = Edata[i] * (1 - Edata[i]*dtime0)/(1 - Edata[i]*dtime)
        if iflag == 3:
            with open('graph_kat.dat', 'w') as f:
                f.write(
                    i, e-18000.0, 
                    (exper-w*(espec[i]-espnm[i]+hnupr*(hnu[i]-espnm2[i]))-backval[i])/error[i],
                    1.0, endef[i]/error[i], w*espnm[i]/error[i], w * espnm2[i]*hnupr/error[i]                
                )
            with open('graph_katd.dat', 'w') as f:
                f.write(
                    i, e-18000.0, exper, w * (espec[i]-espnm[i]+hnupr*(hnu[i]-espnm2[i])) - backval[i],
                    endef[i], w * espnm[i], w * espnm2[i] *hnupr, w * eml[i] * prob_ml
                )

        if e < (emin - 0.01) or e > emax:
            continue
        dif = (exper - theor)/error[i]
        if abs(dif) > cut:
            continue
        fval += dif ** 2
    return fval

    
def specmlint(e, e0, emax):
    # common--------------------
    global thick # source 
    # --------------------------

    # save----------------------
    global ifl_specmlint, de, slev1, cm1
    # --------------------------

    if ifl_specmlint == 0:
        ifl_specmlint += 1
        print("Missing level spectrum calcuation")        # E-=oints nodes list formation
        for i in range(1, 301+1):
            de[i] = i*0.2-2.0
        # End of e-point formation

        # Spline calculation
        acc = .001
        for i in range(301):
            Elow = emax - de[i]
            if Elow >= emax:
                slev1[i] = 0
            else:
                snm = .0
                esp, espnm = specint(Elow, emax, snm, thick)
                accur = esp * acc
                slev1[i] = expmlspectrum(Elow, emax, accur) # ok
        
        spline = interpolate.interp1d(de, slev1, kind = 'cubic') # spline(301, de, slev1, cm1, w1, w2, w3)

        print("Spline is calculated for missing level")
        with open('spn_MCm1.dat', 'w') as f:
            for i in range(301):
                print(de[i], slev1[i], cm1[i])
        
        if e0 - e < 0:
            eml = .0
        else:
            eml = spline(e0-e) # splint(301, e0-e, eml, de, slev1, cm1)
    return eml

def expmlspectrum(e, e0, ac):
    # common-------------
    global pe, pe0 # energy_ml
    # -------------------

    pe = e
    pe0 = e0
    ifail = 0
    ERabs = 1.e-4
    specml = integrate.quad(convol_ml, e, e0, epsabs=ERabs, epsrel=ac) # do1ajf(convol_ml, e, e0, ERabs, ac, specml, er, w, 2000, iw, 260, ifail)
    return specml

def convol_ml(x):
    tspec = truspectrum_ml(x, pe0)
    tran = transmission(x-pe)
    return tspec*tran

def truspectrum_ml(e, tspec, e0):
    pe = np.sqrt(e*(e+2.0))
    de = e0 - e
    de2 = de ** 2
    if de <= 0.:
        tspec = 0.
    else:
        tspec = de2
    tspec = fermi(e) * (e+em)/em*pe*tspec/4.e9
    return tspec

def specint(e, e0, snm, thick):
    # common------------------
    global elev, prolev, ist # levdat
    global shift, shift2 # resfun
    global epoint, npoint # lattice
    # ------------------------

    # save--------------------
    global ifl_specint, de, thickold, shiftold, shift2old, e0old
    global vsnm, vneut1, cm201, cm221, cm021
    global elev1, prolev1, enode, slev, cm1
    # ------------------------


    
    elev1, prolev1 = np.zeros(193), np.zeros(193)
    de, enode, epoint = np.zeros(151), np.zeros(151), np.zeros(301)
    slev1, cm1 = np.zeros(151), np.zeros(151)
    vsnm = np.zeros(110)
    vneut1 = np.zeros(151*110).reshape((151, 110))
    cm201 = np.zeros(151*110).reshape((151, 110))
    cm221 = np.zeros(151*110).reshape((151, 110))
    cm021 = np.zeros(151*110).reshape((151, 110))

    name = ["spn_KNM1.dat"]




    if ifl_specint == 0:
        # read exitation level data
        with open('excitat2.dat', 'r') as f:
            for i in range(201):
                pass

        #-----------
    if abs(shift2 - shift2old) >= .000001:
        print("Final states spectrum with shift=%f" % (shift2))
        slev = 0.
        sl1 = .57408
        sl2 = .42077
        with open("FS_out.dat", "w") as f:
            for i in range(ist):
                if elev1[i] <= 5:
                    prolev[i] = prolev[i] * (sl1 - shift2/100.)/sl1
                else:
                    prolev[i] = prolev[i] * (sl1 - shift2/100.)/sl2
                slev = slev + prolev[i]
                print(i, elev[i], prolev[i], slev, file=f)
        
    ifl_specint += 1
    if (abs(thick-thickold)>=0.000001 or 
        abs(shift2-shift2old)>=0.000001 or 
        abs(shift-shiftold)>=0.000001 or
        abs(e0-e0old)>15.35):

        e0old = e0
        thickold = thick
        shiftold = shift 
        shift2old = shift2
        # calculation of new splines
        print("Needs a spline")
        index = 0
        f = open(name[index], "r")
        read = f.read()
        e0file, thickfile, shiftfile, shift2file = read[0], read[1], read[2], read[3]
        if (abs(thick-thickfile)<=0.0001 and abs(e0-e0file) <= 15.25 and
        abs(shift2-shift2file)<=.0001 and abs(shift-shiftfile)<=.0001):
            print("Reading file=%s with paramters" % (name[i]))
            print("Spectrum endpoint=%f" % (e0file))
            print("Thickness factor= %f" % (thickfile))
            print("Ex/ion shift, %= %f" % (shiftfile))
            print("First state shift, %= %f" % (shift2file))
            for j in range(151):
                de[j], slev1[j], cm1[j] = read[3*j+3], read[3*j+4], read[3*j+5]
            for j in range(110):
                vsnm[j] = read[j+458]
            for j in range(110):
                for k in range(151):
                    vneut1[k][j] = read[568 + 151*j + k]
            f.close()
            
            spline2 = interpolate.interp2d(vneut1, de, vsnm, kind = 'cubic') # spline2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021)
            espnm = spline2(e0-e, snm) # splint2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021, e0-e, snm, espnm, dx, dy, 0, ifail)
            
            spline = interpolate.interp1d(de, slev1, kind = 'cubic') # spline(151, de, slev1, cm1, w1, w2, w3)
            esp = spline(e0-e) # splint(151, e0-e, esp, de, slev1, cm1)

            esp += espnm
            return esp, espnm
        else:
            f.close()
            print("New spline should be calculated")
            print(e0)
            # e-points nodes list formation
            for i in range(9):
                enode[i] = 19261.0-100.0*(i+1)
                print(i, enode[i])
            for i in range(10, 36+1):
                enode[i] = epoint[36-i] + 0.5
                print(i, enode[i])
            for i in range(37, 47+1):
                enode[i] = epoint[0] - 2*(i-37+1)
                print(i, enode[i])
            for i in range(48, 149+1):
                enode[i] = enode[47]-100*(i-48+1)
                print(i, enode[i])
            # end of e-point formation

            # snm point formation
            vsnm[35] = 0.
            for i in range(1, 20+1):
                vsnm[35+i]=0.05*i
                vsnm[35-i]=-vsnm[35+i]
            for i in range(21, 35+1):
                vsnm[35+i]=0.2*(i-20)+1
                vsnm[35-i]=-vsnm[35+i]
            for i in range(71, 109+1):
                vsnm[i]= 4+4*(i-70)**2
            # end of snm point formation

            ac=.0001

            for i in range(151):
                de[i] = e0 - enode[i]
                elow = enode[i]
                vneut1[i][35]=0
                if elow >= e0 - elev1[0]:
                    slev1[i] = 0
                else:
                    snmt = 0.
                    expspectrum(elow, slev1[i], e0, snmt, ac, 3)
                for j in range(35):
                    snmt = vsnm[35+j+1]
                    if elow >= e0 - elev1[0]:
                        snmspec = 0
                    else:
                        expspectrum(elow, snmspec, e0, snmt, 3*ac*slev1[i],7)
                    vneut1[i][35+j+1] = snmspec
                    vneut1[i][35-j-1] = -snmspec
                for j in range(71, 109+1):
                    snmt = vsnm[j]
                    if elow >= e0-elev1[0]:
                        snmspec = 0
                    else:
                        snmspec = expspectrum(elow, snmspec, e0, snmt, 3*ac*slev1[i],7)
                    vneut1[i][j]=snmspec
            
            spline = interpolate.interp1d(de, slev1, kind = 'cubic') # spline(151, de, slev1, cm1, w1, w2, w3)
            spline2 = interpolate.interp2d(vneut1, de, vsnm, kind = 'cubic') # spline2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021)

            print("New splines are calculated with:")
            print("Endpoint energy=%f" % (e0))
            print("Thickness factor= %f" % (thick))
            print("Ex/ion shift, %= %f" % (shift))
            print("First state shift, %= %f" % (shift2))
            with open("spn_last.dat", "w") as f:
                print(e0, thick, shift, shift2)
                for j in range(151):
                    print(de[j], slev1[j], cm1[j], file = f)
                for j in range(110):
                    for k in range(151):
                        print(vneut1[k][j])
            
            espnm = spline2(e0-e, snm) # splint2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021, e0-e, snm, espnm, dx, dy, 0, ifail)
            esp = spline(e0-e) # splint(151, e0-e, esp, de, slev1, cm1)
            esp += espnm
            return esp, espnm


def endeffect(e, eend, step):
    energy = eend - e
    tran = transmission(energy)
    return tran * step

def background(e, e0, back, backpar):
    backval = back + backpar * (18575.0 - e)/40.0
    return

def expspectrum(e, espec, e0, snm, ac, ntyp):
    # common -------------------
    global pe, pe0, psnm, ntype
    # --------------------------

    pe = e
    pe0 = e0
    psnm = snm
    ntype = ntyp
    ifail = 0
    if ntyp < 7:
        erabs = 1.e-10
        espec = integrate.quad(convol, e, e0, epsabs=erabs, epsrel=ac) # d01ajf(convol, e, e0,erabs, ac, espec, er_t, w, 2000, iw, 260, ifail)
    else:
        errel = 1.e-10
        espec = integrate.quad(convol, e, e0, epsabs=ac, epsrel=errel) # d01ajf(convol, e, e0, ac, errel, espec, er_t, w, 2000, iw, 260, ifail)
    return espec

def convol(x):
    # common-------------------------------
    global pe, pe0, psnm, ntype # energy
    # -------------------------------------
    tspec = truspectrum(x, pe0, psnm, ntype)
    tran = transmission(x-pe)
    return tspec*tran

def trapbackground():
    pass

def conlin(x):
    trap = trapsp(x, pe0)
    tran = transmission(x-pe)
    return trap * tran

def trapsp(e, e0):
    pass

def ftrap(x):
    # common ------------------
    global e1, e2 # trapco
    global tspec, scprob # temp
    # -------------------------

    snmt = 0.
    nt = 6
    tspec = truspectrum(x, e2, snmt, nt)
    scprob = scatprob(e1, x, e2)
    return tspec * scprob

def scatprob(e, ep, e0):

    # save----------------
    global ifl_scatprob
    global en, s1, g1
    # --------------------


    ifl_scatprob += 1
    if ifl_scatprob == 1:
        with open("Nscat.dat", "r") as f:
            for i in range(221):
                with open('Nscat.dat', 'r') as f:
                    input(en[i], s1[i])
            spline = interpolate.interp1d(en, s1, kind='cubic') # spline(221, en, s1, g1, w1, w2, w3)
    
    if ep - e < 0.:
        prob = .0
    else:
        prob = spline(ep-e) # splint(221, ep-e, prob, en, s1, g1) #!!!!!
    return prob

def truspectrum(e, e0, snm, ntype):
    # common ----------------------------------
    global elev, prolev, ist # levdat
    # -----------------------------------------

    pe = np.sqrt(e*(e+2*em))
    tspec = 0

    for i in range(ist):
        if ntype == 1: # mainz aprroximation
            de = e0 - elev[i] - e
            de2 = de**2
            if snm >= -1.e-8:
                if de < .0 or de2 < snm:
                    continue
                fuck = de
            if snm < -1.e-8:
                vmu = .76*np.sqrt(-snm)
                if de < -vmu:
                    continue
            fuck = de + vmu*np.exp(-de/vmu - 1)
            tspec = tspec + prolev[i]*fuck*np.sqrt(de2-snm)

        if ntype == 2:
            de = e0 - elev[i] - e
            de2 = de**2
            if de <= .0:
                continue
            tspec = tspec + prolev[i]*(de2 - snm/2)

        if ntype == 3:
            de = e0 - elev[i] - e
            de2 = de**2
            if de <= .0:
                continue
            tspec = tspec + prolev[i]*de2

        if ntype == 6:
            de - e0 - elev[i] - e
            de2 = de**2
            if de <= .0:
                continue
            if de2 > abs(snm):
                dtspec = de*np.sqrt(de2-abs(snm)) - de2
            else:
                dtspec = -de2
            if snm >= .0:
                tspec = tspec + prolev[i]*(de2+dtspec)
            else:
                tspec = tspec + prolev[i]*(de2-dtspec)
        
        if ntype == 7:
            de = e0 - elev[i] - e
            de2 = de**2
            if de <= .0:
                continue
            if de2 > abs(snm):
                dtspec = de*np.sqrt(de2-abs(snm)) - de2
            else:
                dtspec = -de2
            if snm >= 0:
                tspec = tspec + prolev[i]*dtspec
            else:
                tspec = tspec + prolev[i]*(-dtspec)

    if ist == 193 and ntype != 7:
        if e0-e > 165:
            # tail = fstail(e,e0)
            tspec = tspec # + tail
    tspec = fermi(e) * (e+em)/em*pe*tspec/4.e9
    return tspec



def transmission(energy):
    global spline_t
    global shiftold_t, thickold_t

    etrans, trans = np.zeros(221), np.zeros(221)

    if abs(thick - thickold_t) >= 1.0e-6 or abs(shift - shiftold_t >= 1.0e-6):
        thickold_t = thick
        shiftold_t = shift
        RFCAL_KATRIN(thick, etrans, trans)
        spline_t = interpolate.interp1d(etrans, trans, kind='cubic')
    if energy <= 0.:
        return 0.
    if energy >= 300.:
        return 1.
    tran = spline_t(energy) # splint(221, energy, tran, etrans, trans, gtrans)
    return tran

def RFCAL_KATRIN(fac, e, rf, gmr):


    er = np.zeros(221)
    sProb = np.zeros(6)

    sx = 1.
    Bs = 2.52
    Bpinch = 4.24
    sxt = sx * fac
    cos_theta_max = np.sqrt(1 - Bs/Bpinch)
    
    with open('temp.dat', 'w') as f:
        for i in range(100):
            theta_tmp = 1. - .37*np.random.uniform()
            sProb = ProbCalc_int(sxt, sProb, theta_tmp)
            print(theta_tmp, sProb, file=f)

    LOSPAT(er, 0., s2, s2, s3, s4)

    for i in range(221):
        e[i] = er[i]

    erabs = 1.e-5
    rf[0] = 0.
    for i in range(221):
        ifail = 0
        e1 = 0.
        e0 = er[i]
        errel = 1.e-3
        val = integrate.quad(clr, e1, e0, epsabs=erabs, epsrel=errel)
        Spres_point = spres(er[i])
        ProbCalc_int(sxt,sProb,theta_max)
        rf[i]=Spres_point*sProb[0]+val
     
    spline = interpolate.interp1d(er,rf, kind='cubic')
    with open('rf.dat', 'w'):
        for i in range(221):
            print(er[i], rf[i], gmr[i])
    
    return

def ProbCalc_int(x,sp,theta_max):
    if ifl_pc == 0:
        ifl_pc += 1

        P0[39] = 0.
        P1[39] = 0.
        P2[39] = 0.
        P3[39] = 0.
        P4[39] = 0.
        P5[39] = 0.
        theta[39] = 1.

        for i in range(38, -1, -1):
            cos_th = 1 - (39-i)/100 + .005
            theta[i] = 1 - (39-i)/100
            sx = x/cos_th
            spr = ProbCalc(sx, 0, spr)
            P0[i] = P0[i+1] + spr
            spr = ProbCalc(sx, 1, spr)
            P1[i] = P1[i+1] + spr
            spr = ProbCalc(sx, 2, spr)
            P2[i] = P2[i+1] + spr
            spr = ProbCalc(sx, 3, spr)
            P3[i] = P3[i+1] + spr
            spr = ProbCalc(sx, 4, spr)
            P4[i] = P4[i+1] + spr
            spr = ProbCalc(sx, 5, spr)
            P5[i] = P5[i+1] + spr
        
        theta[39] = 1.
        sx = x/theta[39]

        spr = ProbCalc(sx, 0, spr)
        P0[39] = spr
        spr = ProbCalc(sx, 1, spr)
        P1[39] = spr
        spr = ProbCalc(sx, 2, spr)
        P2[39] = spr
        spr = ProbCalc(sx, 3, spr)
        P3[39] = spr
        spr = ProbCalc(sx, 4, spr)
        P4[39] = spr
        spr = ProbCalc(sx, 5, spr)
        P5[39] = spr

        for i in range(40-1):
            P0[i] /= 40 - i
            P1[i] /= 40 - i
            P2[i] /= 40 - i
            P3[i] /= 40 - i
            P4[i] /= 40 - i
            P5[i] /= 40 - i
        
        spline0 = interpolate.interp1d(theta, P0, kind='cubic')
        spline1 = interpolate.interp1d(theta, P1, kind='cubic')
        spline2 = interpolate.interp1d(theta, P2, kind='cubic')
        spline3 = interpolate.interp1d(theta, P3, kind='cubic')
        spline4 = interpolate.interp1d(theta, P4, kind='cubic')
        spline5 = interpolate.interp1d(theta, P5, kind='cubic')

        sp[0] = spline0(theta_max)
        sp[1] = spline1(theta_max)
        sp[2] = spline2(theta_max)
        sp[3] = spline3(theta_max)
        sp[4] = spline4(theta_max)
        sp[5] = spline5(theta_max)

        return





def ProbCalc(x,n,sp):
    spr, pr = np.zeros(6), np.zeros(6)
    pr[0] = np.exp(-x)
    spr[0] = (1. - pr[0])/x

    for i in range(1,n):
        pr[i] = pr[i-1]*x/i
        spr[i] = spr[i-1] - pr[i-1]/i
    
    sp = spr[n-1]
    return

def clr(x):
    global E0, sxt
    sProb = np.zeros(6)
    ep = np.zeros(221)
    Spres_point = spres(E0 - x)
    ProbCalc_int(sxt, sProb, theta_max)
    s1, s2, s3, s4, s5 = LOSPAT(ep,x,s1,s2,s3,s4,s5)
    POG = sProb[0]*s1 + sProb[1]*s2 + sProb[2]*s3 + sProb[3]*s4 + sProb[4]*s5
    return POG * spres(E0 -x)

def spres(E):
    global theta_max

    dE = 0.93
    Bs = 6.
    Bm = 10.
    dE *= 1+shift/100
    if E <= 0:
        spres = 0
        theta_max = 1.
    if E > 0 and E <= dE:
        spres = (1-np.sqrt(1-E/dE*Bs/Bm))/(1-np.sqrt(1-Bs/Bm))
        theta_max = np.sqrt(1-E/dE*Bs/Bm)
    if E > dE:
        spres = 1
        theta_max = np.sqrt(1-Bs/Bm)
    return spres


def LOSPAT(ep, e_point, s1p, s2p, s3p, s4p, s5p):

    global e, s1, s2, s3, s4, s5

    if ifl_ls == 0:
        ifl_ls += 1

        for i in range(21):
            e[i] = (i-1)/20.0
        
        for i in range(21, 190+1):
            e[i] = .2*(i-20)+1

        for i in range(191, 211+1):
            e[i] = 3. *(i-190) + 35

        for i in range(212, 220+1):
            e[i] = 20.*(i-211) + 100
        


        xc1 = 12.6
        xc2 = 14.3
        w_1 = 1.85
        w_2 = 12.5
        a1 = .204
        a2 = .0556

        for i in range(221):
            de = e[i]
            if de <= 14.3:
                s1[i] = a1 * np.exp(-2.*(de-xc1)**2/w_1**2)
            else:
                s1[i] = a2 * w_2**2/(w_2**2 + 4.*(de-xc2)**2)
        
        spline1 = interpolate.interp1d(e, s1, kind='cubic')

        ifail = 0
        
        # 2-scatterings
        s2[0] = 0.
        for i in range(1, 220+1):
            Emax = e[i]
            erabs = 1.e-10
            errel = 1.e-4
            Emin = 0.
            C11 = lambda x: spline1(x) * spline1(Emax - x)
            s2[i] = integrate.quad(C11, Emin, Emax, epsabs = erabs, epsrel=errel)
        
        spline2 = interpolate.interp1d(e, s2, kind='cubic')

        # 3-scatterings
        s3[0] = 0.
        for i in range(1, 220+1):
            Emax = e[i]
            C12 = lambda x: spline1(x) * spline2(Emax - x)
            s3[i] = integrate.quad(C12, Emin, Emax, epsabs=erabs, epsrel=errel)
        
        spline3 = interpolate.interp1d(e, s3, kind='cubic')

        # 4-scatterings
        s4[i] = 0.
        for i in range(1, 220+1):
            Emax = e[i]
            C13 = lambda x: spline1(x) * spline3(Emax - x)
            s4[i] = integrate.quad(C13, Emin, Emax, epsabs=erabs, epsrel=errel)

        spline4 = interpolate.interp1d(e, s4, kind='cubic')

        # 5-scatterings
        s5[i] = 0.
        for i in range(1, 220+1):
            Emax = e[i]
            C14 = lambda x: spline1(x) * spline4(Emax - x)
            s5[i] = integrate.quad(C14, Emin, Emax, epsabs=erabs, epsrel=errel)

        spline5 = interpolate.interp1d(e, s5, kind='cubic')

    s1p = spline1(e_point)
    s2p = spline2(e_point)
    s3p = spline3(e_point)
    s4p = spline4(e_point)
    s5p = spline5(e_point)

    return s1p, s2p, s3p, s4p, s5p



def fermi(e):
    beta = np.sqrt(1-(em/(em+e))**2)
    y = 4*np.pi*alpha/beta
    return y/abs(1 - np.exp(-y)) * (1.002037 - 0.001427 * beta)
