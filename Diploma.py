import numpy as np
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
slev1, cm1 = np.zeros(151), np.seros(151)
vsnm = np.zeros(110)
vneut1 = np.zeros(151*110).reshape((151, 110))
cm201 = np.zeros(151*110).reshape((151, 110))
cm221 = np.zeros(151*110).reshape((151, 110))
cm021 = np.zeros(151*110).reshape((151, 110))

# scatprob------------------------
en, s1, g1 = np.zeros(221), np.zeros(221), np.zeros(221)
ifl_scatprob = 0

#------------------------------------------------





def chi2(npar, grad, fval, xval, iflag):
    # common------------------
    global thick # source
    global shift, shift2 # resfun
    global epoint, npoint # lattice
    # ------------------------

    # save--------------------
    global ifl_chi2
    global Edata, error
    # ------------------------


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
                i+=1
                print("%3d %10.1f %11.6f %13.8f\n" % (i, epoint[i], Edata[i], error[i]))

        
        # Arrays calculation
        for i in range(1, 101):
            entr = float(i)
            tran = transmission(entr)
            # fstail(e0 - entr - 165., e0, tail)
            specint(e0 - entr, vesp, e0, snm, vespnm, thick)
            print(entr, tran)
        # Arrays calculated

        # Faked spectrum formation--------------------------------------
        w = 96.4
        with open('run_kat_MC2.dat', 'w') as f:
            print("Simulation with Elow=%f" % (emin), file=f)
            print("%s %s %s %s" % ("HV", "Freq", "Err", "Time"), file=f)
            for i in range(npoint):
                e = epoint[i]
                specint(e, espec[i], e0, snm, espnm[i], thick)
                specint(e, hnu[i], e0, snm2, espnm2[i], thick)
                specmlint(e, eml[i], pos_ml, e0)
                trapbackground(e, trspec[i], e0)
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
            specint(e, espec[i], e0, snm, espnm[i], thick)
            specint(e, hnu[i], e0, snm, espnm[i], thick)
        if ntype == 2:
            expspectrum(e, espec[i], e0, snm, acc)
            expspectrum(e, hnu[i], e0, snm2, acc)
        specmlint(e, eml[i], pos_ml, e0)
        background(e, e0, back, backpar, backval[i])
        endef[i] = endeffect(e, eend, step)
        # if
        sespec = sespec + espec[i] + hnu[i]*hnupr
        sml = sml +eml[i]*prob_ml
        sbackval = sbackval + backval[i]
        sendef = sendef + endef[i]
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

    
def specmlint(e, eml, e0, emax):
    # common--------------------
    global thick # source 
    # --------------------------

    # save----------------------
    global ifl_specmlint, de, slev1, cm1
    # --------------------------

    w1, w2, w3 = np.zeros(301), np.zeros(301), np.zeros(301)


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
                specint(Elow, esp, emax, snm, espnm, thick)
                accur = esp * acc
                expmlspectrum(Elow, slev1[i], emax, accur)
        
        spline(301, de, slev1, cm1, w1, w2, w3)

        print("Spline is calculated for missing level")
        with open('spn_MCm1.dat', 'w') as f:
            for i in range(301):
                print(de[i], slev1[i], cm1[i])
        
        if e0 - e < 0:
            eml = .0
        else:
            splint(301, e0-e, eml, de, slev1, cm1)
    return

def expmlspectrum(e, specml, e0, ac):
    # common-------------
    global pe, pe0 # energy_ml
    # -------------------

    w, iw = np.zeros(2000), np.zeros(260)
    pe = e
    pe0 = e0
    ifail = 0
    ERabs = 1.e-4
    do1ajf(convol_ml, e, e0, ERabs, ac, specml, er, w, 2000, iw, 260, ifail)

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

def specint(e, esp, e0, snm, espnm, thick):
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


    
    w1, w2, w3 = np.zeros(151), np.zeros(151), np.zeros(151)
    elev1, prolev1 = np.zeros(193), np.zeros(193)
    de, enode, epoint = np.zeros(151), np.zeros(151), np.zeros(301)
    slev1, cm1 = np.zeros(151), np.seros(151)
    vsnm = np.zeros(110)
    vneut1 = np.zeros(151*110).reshape((151, 110))
    cm201 = np.zeros(151*110).reshape((151, 110))
    cm221 = np.zeros(151*110).reshape((151, 110))
    cm021 = np.zeros(151*110).reshape((151, 110))

    name = ["spn_KNM1.dat"]




    if ifl_specint == 0:
        # read exitation level data

        #-----------
        pass
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
            spline(151, de, slev1, cm1, w1, w2, w3)
            spline2(151, 110, vneut, de, vsnm, cm201, cm221, cm021)
            splint2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021, e0-e, snm, espnm, dx, dy, 0, ifail)
            splint(151, e0-e, esp, de, slev1, cm1)
            esp += espnm
            return
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
                        expspectrum(elow, snmspec, e0, snmt, 3*ac*slev1[i],7)
                    vneut1[i][j]=snmspec
            
            spline(151, de, slev1, cm1, w1, w2, w3)
            spline2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021)

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
            
            splint2(151, 110, vneut1, de, vsnm, cm201, cm221, cm021, e0-e, snm, espnm, dx, dy, 0, ifail)
            splint(151, e0-e, esp, de, slev1, cm1)
            esp += espnm
            return


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
        d01ajf(convol, e, e0,erabs, ac, espec, er_t, w, 2000, iw, 260, ifail)
    else:
        errel = 1.e-10
        d01ajf(convol, e, e0, ac, errel, espec, er_t, w, 2000, iw, 260, ifail)
    return

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
    conlin = trap * tran

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

    w1, w2, w3 = np.zeros(221), np.zeros(221), np.zeros(221)

    ifl_scatprob += 1
    if ifl_scatprob == 1:
        with open("Nscat.dat", "r") as f:
            for i in range(221):
                pass
            spline(221, en, s1, g1, w1, w2, w3)
    
    if ep - e < 0.:
        prob = .0
    else:
        splint(221, ep-e, prob, en, s2, g1) #!!!!!
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
            tail = fstail(e,e0)
            tspec = tspec + tail
    tspec = fermi(e) * (e+em)/em*pe*tspec/4.e9
    return tspec



def transmission(energy):
    pass

def RFCAL_KATRIN(fac, e, rf, gmr):
    pass

def fermi(e):
    beta = np.sqrt(1-(em/(em+e))**2)
    y = 4*np.pi*alpha/beta
    return y/abs(1 - np.exp(-y)) * (1.002037 - 0.001427 * beta)
