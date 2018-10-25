import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pprint
#matplotlib.use('Agg')
# copy number
c = 5
# cooperativity of repressor binding
n = 1
# transcription rates
aR = 1.03
aGmax = 3.78
aGmin = 0.0001
# degradation rates
yR = 0.06 #estimated
yG = 0.018
yMR = 0.2
yMG = 0.2
# translation rates
bR = 0.54
bG = 3.647
# on and off rates of repressor binding to the promoter
kRon = 50
kRoff = 334
kDn = kRoff / kRon

def wrap(c, aR, aGmax, aGmin, yR, yG, yMR, yMG, bR, bG, kRon, kRoff):
    def model(z, t):
        mR = z[0]
        R = z[1]
        PG = z[2]
        PGR = z[3]
        mG = z[4]
        G = z[5]
        Rn=R
        dmRdt = c * aR - yMR * mR
        # dRdt = bR * mR - yR * R - n * kRon * Rn * PG + n * kRoff * PGR + (n - 1) * n * yR * PGR
        dRdt = bR * mR - yR * R #+ kRoff * PGR - kRon * Rn * PG
        dPGdt = kRoff * PGR - kRon * Rn**n * PG + n * yR * PGR
        dPGRdt = kRon * Rn**n * PG - kRoff * PGR - n * yR * PGR
        dmGdt = aGmax * PG + aGmin * PGR - yMG * mG
        # dmGdt = c * (aGmin + (aGmax - aGmin) * (kDn / (kDn + Rn))) - yM * mG
        # dmGdt = c * (aGmax * kDn / Rn) - yM * mG
        dGdt = bG * mG - yG * G
        dzdt = [dmRdt, dRdt, dPGdt, dPGRdt, dmGdt, dGdt]
        return dzdt

    z0 = [0, 0, c, 0, 0, 0]
    t = np.linspace(0, 1000, 10000)
    return odeint(model, z0, t)

#def main(c, par, kRon, kRoff, aGmax): #kRon, kRoff, yR,yG,yMR, yMG, aR,aGmax,bR,bG
def main(c,par):
    goi=[]
    for c in range(c[0],c[1]+1,c[2]):
        z = wrap(c, aR, aGmax, aGmin, yR, yG, yMR, yMG, bR, bG, kRon,kRoff)
        goi.append(z[:,5][-1])
        if True:
            t = np.linspace(0, 1000, 10000)
            plt.plot(t,z[:,0],'b-',label='mR')
            plt.plot(t,z[:,1],'r-',label='R')
            plt.plot(t,z[:,2],'b-.',label='PG')
            plt.plot(t,z[:,3],'r-.',label='PGR')
            plt.plot(t,z[:,4],'b--',label='mG')
            plt.plot(t,z[:,5],'g-',label='G')
            plt.ylabel('concentration')
            plt.xlabel('time')
            plt.legend(loc='best')
            plt.show()
    return goi



if __name__ == '__main__':
    r = main([5,100,5],'df')
    print(r)
