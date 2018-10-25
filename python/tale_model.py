import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import pprint

# copy number
c = 30
# cooperativity of repressor binding
n = 1
# transcription rates
aR = 1.03
aGmax = 3.78
aGmin = 0.0001
# degradation rates
yR = 0.047 #estimated
yG = 0.019
yM = 0.2
# translation rates
bR = 0.55
bG = 3.647
# on and off rates of repressor binding to the promoter
kRon = 1.6
kRoff = 3.85
kDn = kRoff / kRon

def wrap(c, aR, aGmax, aGmin, yR, yG, yM, bR, bG):
    def model(z, t):
        mR = z[0]
        R = z[1]
        PG = z[2]
        PGR = z[3]
        mG = z[4]
        G = z[5]

        #print(aGmax*kDn/(kDn+R))
        #input()

        Rn=R
        dmRdt = c * aR - yM * mR
        # dRdt = bR * mR - yR * R - n * kRon * Rn * PG + n * kRoff * PGR + (n - 1) * n * yR * PGR
        dRdt = bR * mR - yR * R #+ kRoff * PGR - kRon * Rn * PG
        dPGdt = kRoff * PGR - kRon * Rn**n * PG + n * yR * PGR
        dPGRdt = kRon * Rn**n * PG - kRoff * PGR - n * yR * PGR
        dmGdt = aGmax * PG + aGmin * PGR - yM * mG
        # dmGdt = c * (aGmin + (aGmax - aGmin) * (kDn / (kDn + Rn))) - yM * mG
        # dmGdt = c * (aGmax * kDn / Rn) - yM * mG
        dGdt = bG * mG - yG * G
        dzdt = [dmRdt, dRdt, dPGdt, dPGRdt, dmGdt, dGdt]
        return dzdt

    z0 = [0, 0, c, 0, 0, 0]
    t = np.linspace(0, 1000, 10000)
    return odeint(model, z0, t)

c_number_min = 5
c_number_max = 101
c_number_step = 5

def main():
    final_results = []
    final_R=[]
    final=[]
    for c in range(c_number_min, c_number_max, c_number_step):
        z = wrap(c, aR, aGmax, aGmin, yR, yG, yM, bR, bG)
        final_results.append(z[:,5][-1])
        final_R.append(z[:,1][-1])
        if False:
            t = np.linspace(0, 1000, 10000)
            # plt.plot(t,z[:,0],'b-',label='mR')
            #plt.plot(t,z[:,1],'r-',label='R')
            #plt.plot(t,z[:,2],'b-.',label='PG')
            #plt.plot(t,z[:,3],'r-.',label='PGR')
            # plt.plot(t,z[:,4],'b--',label='mG')
            plt.plot(t,z[:,5],'g-',label='G')
            #plt.plot( z[:,1],aGmax*(kDn/(kDn+z[:,1])),'b--')
            #plt.ylabel('concentration')
            #plt.xlabel('time')
            #plt.legend(loc='best')
            plt.show()
    if True:
        error=(final_results[-1] - final_results[0]) / final_results[0]
        ax=plt.axes()
        line=plt.plot(range(c_number_min, c_number_max, c_number_step),final_results,'xkcd:cornflower blue')

        #plt.yscale("log")
        #plt.xscale("log")
        #plt.ylim(400,1000)
        #ax.set_yticks([10,100,1000])
        #ax.tick_params(labelleft=True,bottom=True,left=True,top=True,right=True, direction='in')
        #plt.text(10,500,'Error= ' + str(error))
        #print(ax.get_yaxis())
        #ax.get_yaxis().set_major_formatter(LogFormatterMathtext())
        #plt.plot(range(c_number_min, c_number_max, c_number_step),range(c_number_min, c_number_max, c_number_step),'r-.',label='Reality')
        plt.ylabel('GOI')
        plt.xlabel('Copy Number')
        #plt.legend('E: ' + str(error),loc='center')
        plt.show()
        #plt.plot(final_R,final,'r-',label='GOI-R')
        #plt.ylabel('rate')
        #plt.xlabel('R')
        #plt.show()

    #beta = bR * aR / (yR * yM)
    #error = kDn / beta * c_number_min
    print((final_results[-1] - final_results[0]) / final_results[0])
    #pprint.pprint(final_results)
    #strength = kDn / (kDn + beta * c_number_min)
    #return {
    #    'beta': beta,
    #    'error': error,
    #    'strength': strength,
    #    'min_goi': min(final_results),
    #    'max_goi': max(final_results),
    #}
    write(range(c_number_min, c_number_max, c_number_step),final_results)

def write(c,goi):
    name= "goi_c_TALE"
    filename = name + '.csv'
    columns='c,goi'
    file = open(filename, 'w')
    file.write(columns + '\n')
    for i in range(len(c)):
        file.write(str(c[i]) + ',' + str(goi[i]) + '\n')

if __name__ == '__main__':
    main()
