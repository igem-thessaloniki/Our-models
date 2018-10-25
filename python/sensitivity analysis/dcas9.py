import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# copy number
c = 5
cn = 100
genome = 1
n = 2

Dox=24.4
Rham=100000

at=4.33
ah=3.22

bt=1.819
bd=0.55
bh=1.304
bG=3.65

dmt=0.2
dt=0.017
dmd=0.2
dd=0.005 #(0.001-1)
dmh=0.2
dh=0.06 #(0.001-1)
ds=0.18 #(0.001-1)
dmG=0.2
dG=0.0193

Tmax=0.73
Tmin=0.0001
hmax=5
hmin=0.0001
aGmax=3.78
aGmin=0.0001

kdtoff=0.12#(0.001-10)
kdton=126#(0.001-10)
kdson=0.11#(0.001-10)
kdsoff=0.001#(0.001-10)
khron=100#(0.001-10)
khroff=100#(0.001-10)
kon1=110#(0.001-10)
koff1=1#(0.001-10)
kon2=0.000001#(0.001-10)
koff2=100#(0.001-10)
kon3=0.693#(0.001-10)
koff3=0.18#(0.001-10)

def wrap(par, kon3,koff3,kdson,kdsoff,kdton,kdtoff,khron,khroff,kon1,koff1, kon2, koff2, Dox, at, bt, bd,
        dmt, dt, dmd, dd, Tmax, Rham,ah,bh,dmh,dh,ds,hmax):
    def model(z, t):
        mRtetR=z[0]
        tetR=z[1]
        Dox=z[2]
        Dox_tetR=z[3]
        P=z[4]
        P_tetR=z[5]
        mRcas9 = z[6]
        cas9 = z[7]


        mRrhaS=z[8]
        rhaS=z[9]
        Rham = z[10]
        rhaS_Rham = z[11]
        Pbad = z[12]
        Pbad_rhaS_Rham = z[13]
        sgRNA = z[14]

        R = z[15]
        PG = z[16]
        PG_R = z[17]
        mG = z[18]
        G = z[19]

        dmRtetRdt = at * cn - dmt * mRtetR
        dtetRdt = bt * mRtetR - dt * tetR + kdtoff * Dox_tetR - kdton * Dox * tetR + koff2 * P_tetR - kon2 * P * tetR**n

        dDoxdt =  kdtoff * Dox_tetR - kdton * Dox * tetR + dt * Dox_tetR
        if Dox<0: Dox=0
        if Dox>24.4: Dox=24.4
        dDox_tetRdt = kdton * Dox * tetR - kdtoff * Dox_tetR - dt * Dox_tetR

        dPdt = koff2 * P_tetR - kon2 * P * (tetR ** n) + n* dt * P_tetR
        dP_tetRdt = kon2 * P * (tetR**n) - koff2 * P_tetR - n * dt * P_tetR
        dmRcas9dt = Tmax * P + Tmin * P_tetR - dmd * mRcas9
        if par == 'dcas' or par == 'dcasred': dcas9dt = bd * mRcas9 - dd * cas9 #+ kdsoff * R - kdson * sgRNA * cas9
        else:dcas9dt = bd * mRcas9 - dd * cas9 + kdsoff * R - kdson * sgRNA * cas9

        dRhamdt =  -khron * rhaS * Rham + khroff * rhaS_Rham + dh * rhaS_Rham
        dmRrhaSdt = ah * genome - dmh * mRrhaS
        drhaSdt = bh * mRrhaS - dh * rhaS - khron * rhaS * Rham + khroff * rhaS_Rham
        if Rham<0 : Rham=0
        if Rham>1210000: Rham=1210000
        drhaS_Rhamdt = khron * rhaS * Rham - khroff * rhaS_Rham - dh * rhaS_Rham - kon1 * rhaS_Rham * Pbad + koff1 * Pbad_rhaS_Rham + dh * Pbad_rhaS_Rham
        dPbaddt = -kon1 * rhaS_Rham * Pbad + koff1 * Pbad_rhaS_Rham + dh * Pbad_rhaS_Rham
        dPbad_rhaS_Rhamdt = -koff1 * Pbad_rhaS_Rham + kon1 * rhaS_Rham * Pbad - dh * Pbad_rhaS_Rham
        if par=='sgrna' : dsgRNAdt= hmax * Pbad_rhaS_Rham + hmin * Pbad - ds * sgRNA# - kdson * sgRNA * cas9 + kdsoff * R
        else: dsgRNAdt= hmax * Pbad_rhaS_Rham + hmin * Pbad - ds * sgRNA - kdson * sgRNA * cas9 + kdsoff * R

        dRdt = kdson * cas9 * sgRNA - kdsoff * R - dd * R - kon3 * R * PG + koff3 * PG_R
        dPGdt = koff3 * PG_R - kon3 * PG * R + dd * PG_R
        dPG_Rdt = -koff3 * PG_R + kon3 * PG * R - dd * PG_R
        dmGdt = aGmax * PG + aGmin * PG_R - dmG * mG
        dGdt = bG * mG - dG * G
        dzdt = [dmRtetRdt, dtetRdt, dDoxdt, dDox_tetRdt, dPdt, dP_tetRdt, dmRcas9dt, dcas9dt, dmRrhaSdt, drhaSdt, dRhamdt, drhaS_Rhamdt, dPbaddt, dPbad_rhaS_Rhamdt, dsgRNAdt, dRdt, dPGdt, dPG_Rdt, dmGdt, dGdt]
        return dzdt

    z0 = [0,0,Dox,0,temp,0,0,0,0,0,Rham,0,c,0,0,0,c,0,0,0]
    t = np.linspace(0, 2000, 4000)
    return odeint(model, z0, t)

def main(c,par,kdson,kdsoff,kon3,koff3,kon1,Rham):
#def main(par):
    z = wrap(par, kon3,koff3,kdson,kdsoff,kdton,kdtoff,khron,khroff,kon1,koff1, kon2, koff2, Dox, at, bt,bd,
            dmt, dt, dmd, dd, Tmax, Rham,ah,bh,dmh,dh,ds,hmax)
    if False:
        t = np.linspace(0, 2000, 4000)
        #plt.plot(t,z[:,1],'b-.',label='TetR')
        #plt.plot(t,z[:,2],'b--.',label='Dox')
        #plt.plot(t,z[:,3],'r-.',label='Dox_tetR')
        #plt.plot(t,z[:,2],'g-.',label='PGR1')
        #plt.plot(t,z[:,4],'r--.',label='PG2')
        #plt.plot(t,z[:,13],'g--.',label='Pbad2')
        #plt.plot(t,z[:,12],'b-',label='Pbad')
        #plt.plot(t,z[:,7],'r--',label='cas9')
        #plt.plot(t,z[:,14],'g-',label='sgRNA')
        plt.plot(t,z[:,15],'b-.',label='R')
        #plt.plot(t,z[:,10],'r-.',label='Rham')
        #plt.plot(t,z[:,9],'g-.',label='Rhas')
        #plt.plot(t,z[:,11],'b--',label='simploko')
        plt.plot(t,z[:,19],'r--',label='G')
        plt.ylabel('concentration')
        plt.xlabel('time')
        plt.legend(loc='best')
        plt.show()


    if par== 'dcas' or par=='dcasred': return z[:,7][-1]
    elif par == 'sgrna' : return z[:,14][-1]
    else: return z[:,19][-1]



if __name__ == '__main__':
    r = main('dcas')
    print(r)
