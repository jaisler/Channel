import numpy as np
import pandas as pd
import math
import sys
import subprocess

class DataPSD:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = [] # attribute or member
        self.v = []
        self.w = []
        self.umag = []
        self.umag2 = []
        self.p = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfilepsd']):
            #subprocess.call(["sed -e 's/ /,/g' "
            #          + params['path'] + "/"
            #          + params['filePSD'] + str(i)+ ".his > "
            #          + params['path'] + "/"
            #          + params['filePSD'] + str(i) + ".csv"], shell=True)
            df.append(pd.read_csv(params['path'] + '/'
                + params['filePSD'] + str(i) + '.csv',
                delimiter=',',
                skiprows=[j for j in range(1, params['psdloc'])]))

            # Obtain data in each psdNpoints rows
            df[i] = df[i].iloc[::params['psdNpoints'], :]

            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rho'])
            self.umag.append(np.sqrt((df[i]['rhou']/df[i]['rho'])**2 
                                   + (df[i]['rhov']/df[i]['rho'])**2 
                                   + (df[i]['rhow']/df[i]['rho'])**2))
            self.umag2.append(((df[i]['rhou']/df[i]['rho'])**2 
                                   + (df[i]['rhov']/df[i]['rho'])**2 
                                   + (df[i]['rhow']/df[i]['rho'])**2))
            self.p.append((params['gamma'] - 1) * (df[i]['E'] - 0.5 * df[i]['rho'] * self.umag2[i]))
            self.t.append(df[i]['t'])
            self.dt.append(df[i]['t'][params['psdNpoints']]-df[i]['t'][0])
            
            # Calculates the mean velocity and then gets the fluctuation
            tn = 0
            denu = 0
            denv = 0
            denw = 0
            denumag = 0
            denp = 0
            for j in range(len(self.u[i])):
                denu += self.u[i][j*params['psdNpoints']]
                denv += self.v[i][j*params['psdNpoints']]
                denw += self.w[i][j*params['psdNpoints']]
                denumag += self.umag[i][j*params['psdNpoints']]
                denp += self.p[i][j*params['psdNpoints']]
                tn += 1
            denu /= tn
            denv /= tn
            denw /= tn
            denumag /= tn
            denp /= tn

            self.u[i] = self.u[i] - denu
            self.v[i] = self.v[i] - denv
            self.w[i] = self.w[i] - denw
            self.w[i] = self.umag[i] - denumag
            self.w[i] = self.p[i] - denp

    def GetU(self):
        return self.u

    def GetV(self):
        return self.v

    def GetW(self):
        return self.w

    def GetUmag(self):
        return self.umag

    def GetP(self):
        return self.p

    def GetT(self):
        return self.t

    def GetDt(self):
        return self.dt
    
class DataSpatialCorrelation:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = [] # attribute or member
        self.v = []
        self.w = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfilecorr']):
            #subprocess.call(["sed -e 's/ /,/g' "
            #          + params['path'] + "/"
            #          + params['fileCorr'] + str(i)+ ".his > "
            #          + params['path'] + "/"
            #          + params['fileCorr'] + str(i) + ".csv"], shell=True)
            df.append(pd.read_csv(params['path'] + '/'
                + params['fileCorr'] + str(i) + '.csv',
                delimiter=',',))
            
            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rho'])
            self.t.append(df[i]['t'])
            self.dt.append(df[i]['t'][params['psdNpoints']]-df[i]['t'][0])

            # Calculates the mean velocity and then gets the fluctuation
            k = 0
            tn = 0
            denu = 0
            denv = 0
            denw = 0
            while (k < len(self.u[i])):
                meanuz = 0
                meanvz = 0
                meanwz = 0
                for j in range(params['npcorr']): 
                    meanuz += self.u[i][j + k]
                    meanvz += self.v[i][j + k]
                    meanwz += self.w[i][j + k]
                meanuz /= (params['npcorr'])
                meanvz /= (params['npcorr'])
                meanwz /= (params['npcorr'])
                denu += meanuz
                denv += meanvz
                denw += meanwz
                k += params['npcorr']
                tn += 1
            denu /= tn
            denv /= tn
            denw /= tn

            self.u[i] = self.u[i] - denu
            self.v[i] = self.v[i] - denv
            self.w[i] = self.w[i] - denw

    def GetUprime(self):
        return self.u

    def GetVprime(self):
        return self.v

    def GetWprime(self):
        return self.w

    def GetT(self):
        return self.t

    def GetDt(self):
        return self.dt
    
class DataTBLVProfile:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = []
        self.v = []
        self.w = []
        self.p = []
        self.umag2 = []
        self.urmsp = []
        self.vrmsp = []
        self.wrmsp = []
        self.prmsp = []
        self.uvp = []
        self.upp = []
        self.uplus = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfileprof']):
            #subprocess.call(["sed -e 's/ /,/g' "
            #          + params['path'] + "/"
            #          + params['fileProf'] + str(i)+ ".his > "
            #          + params['path'] + "/"
            #          + params['fileProf'] + str(i) + ".csv"], shell=True)
            df.append(pd.read_csv(params['path'] + '/'
                + params['fileProf'] + str(i) + '.csv',
                delimiter=',',))

            # Prepare data: u,v,w
            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rho'])
            self.umag2.append(((df[i]['rhou']/df[i]['rho'])**2 
                                   + (df[i]['rhov']/df[i]['rho'])**2 
                                   + (df[i]['rhow']/df[i]['rho'])**2))
            self.p.append((params['gamma'] - 1) * (df[i]['E'] - 0.5 * df[i]['rho'] * self.umag2[i]))

            # U+
            up = np.zeros(params['npointsy'])
            l = 0
            tn = 0
            while(l < len(self.u[i])):
                for k in range(params['npointsy']):
                    uAvgZ = 0
                    for j in range(params['npointsz']):
                        # it is divided by the u_{tau} to obtain U+
                        uAvgZ += self.u[i][j + k * params['npointsz'] + l]
                    # Obtain the average in z-direction
                    uAvgZ = uAvgZ / params['npointsz']
                    # Obtain the average in time
                    up[k] += uAvgZ 
                l += params['npointsz'] * params['npointsy']            
                tn += 1
            up /= (tn * params['utau']) # U/u_{tau}

            # P mean
            pm = np.zeros(params['npointsy'])
            l = 0
            tn = 0
            while(l < len(self.p[i])):
                for k in range(params['npointsy']):
                    pAvgZ = 0
                    for j in range(params['npointsz']):
                        pAvgZ += self.p[i][j + k * params['npointsz'] + l]
                    # Obtain the average in z-direction
                    pAvgZ = pAvgZ / params['npointsz']
                    # Obtain the average in time
                    pm[k] += pAvgZ 
                l += params['npointsz'] * params['npointsy']            
                tn += 1
            pm /= tn 

            #RMSplus
            qunt = params['qunt']
            rmspz = np.zeros((qunt, params['npointsy']))
            for j in range(params['npointsz']): # z
                rmsp = np.zeros((qunt, params['npointsy']))
                tn = 0
                l = 0
                while(l < len(self.u[i])): # time
                    for k in range(params['npointsy']): # y
                        rmsp[0][k] += (self.u[i][j + k * params['npointsz'] + l] - up[k] * params['utau'])**2
                        rmsp[1][k] += self.v[i][j + k * params['npointsz'] + l]**2
                        rmsp[2][k] += self.w[i][j + k * params['npointsz'] + l]**2
                        rmsp[3][k] += ((self.u[i][j + k * params['npointsz'] + l] - up[k] * params['utau']) 
                                      * self.v[i][j + k * params['npointsz'] + l])
                        rmsp[4][k] += (self.p[i][j + k * params['npointsz'] + l] - pm[k])**2
                        rmsp[5][k] += ((self.u[i][j + k * params['npointsz'] + l] - up[k] * params['utau']) 
                                      * (self.p[i][j + k * params['npointsz'] + l] - pm[k]))

                    l += params['npointsy'] * params['npointsz']          
                    tn += 1
            
                for l in range(qunt-3):
                    rmsp[l] /= tn 
                    rmsp[l] = np.sqrt(rmsp[l]) / params['utau'] # <uu>^{1/2}/u_{\tau}
                    rmspz[l] += rmsp[l]                 
                rmsp[qunt-3] /= (tn * (params['utau']**2)) # <uv>/u_{\tau}^2
                rmspz[qunt-3] += rmsp[qunt-3]
                # pressure related quantities
                # p^{+}_{rms}
                rmsp[qunt-2] /= tn 
                rmsp[qunt-2] = np.sqrt(rmsp[qunt-2]) / (params['rhom']*(params['utau'])**2)
                rmspz[qunt-2] += rmsp[qunt-2]
                # <up> / (\rho * u_{tau}^3)
                rmsp[qunt-1] /= tn 
                rmsp[qunt-1] = rmsp[qunt-1] / (params['rhom']*(params['utau'])**3)
                rmspz[qunt-1] += rmsp[qunt-1]

            for i in range(qunt):
                rmspz[i] /= params['npointsz']

            self.uplus.append(up)
            self.urmsp.append(rmspz[0])
            self.vrmsp.append(rmspz[1])
            self.wrmsp.append(rmspz[2])
            self.uvp.append(-rmspz[3]) 
            self.prmsp.append(rmspz[4])
            self.upp.append(rmspz[5])


    def GetUplus(self):
        return self.uplus

    def GetURMSP(self):
        return self.urmsp

    def GetVRMSP(self):
        return self.vrmsp

    def GetWRMSP(self):
        return self.wrmsp   

    def GetUVP(self):
        return self.uvp
    
    def GetPRMSP(self):
        return self.prmsp  

    def GetUPP(self):
        return self.upp 
    
class DataDNSDataReyStress:
    def __init__(self, params):
        """ Extract DNS data: Reynolds Stresses """

        self.uu = [] # attribute or member
        self.vv = []
        self.ww = []
        self.uv = []
        self.y = []
        self.yplus = []
        drey = pd.read_csv(params['pathE'] + '/'
            + params['fileRSE'] + '.csv', delimiter=',')

        self.uu = np.sqrt(drey['uu'])
        self.vv = np.sqrt(drey['vv'])
        self.ww = np.sqrt(drey['ww'])
        self.uv = -drey['uv']
        self.y = drey['y']
        self.yplus = drey['yplus']

    def GetURMS(self):
        return self.uu
    
    def GetVRMS(self):
        return self.vv

    def GetWRMS(self):
        return self.ww

    def GetUV(self):
        return self.uv

    def GetY(self):
        return self.y
    
    def GetYplus(self):
        return self.yplus

class DataDNSDataProfile:
    def __init__(self, params):
        """ Extract DNS data: Profile """

        self.umean = [] # attribute or member
        self.y = []
        self.yplus = []
        dprof = pd.read_csv(params['pathE'] + '/'
            + params['filePrE'] + '.csv', delimiter=',')

        self.umean = dprof['Umean']
        self.y = dprof['y']
        self.yplus = dprof['yplus']

    def GetUmean(self):
        return self.umean

    def GetY(self):
        return self.y

    def GetYplus(self):
        return self.yplus

class DataDNSDataProfilePressure:
    def __init__(self, params):
        """ Extract DNS data: Profile """
        
        self.prmsp = [] # attribute or member
        self.upp = []
        self.yplus = []
    
        #for i in range(params['nfilesP']):
        dp = pd.read_csv(params['pathE'] + '/'
                + params['filePP'] + '.csv', delimiter=',')

        self.prmsp = np.sqrt(dp['R_pp'])
        self.upp = dp['R_up']
        self.yplus = dp['yplus']

    def GetPRMSP(self):
        return self.prmsp

    def GetUPP(self):
        return self.upp

    def GetYplus(self):
        return self.yplus