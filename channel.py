import numpy as np
import math
import yaml

from classes import *
import functions as fc
import plot as pl

with open(r'configuration.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params = yaml.load(file, Loader=yaml.FullLoader)

if params['post']:
    if params['routine']['psd']:
        # Initialisation: data from the probes.
        objPSD = DataPSD(params)
        dt = objPSD.GetDt()
        # Pressure
        p = objPSD.GetP()
        (fp, Sp) = fc.CalculatePSD(p, dt, params)
        for i in range(len(fp)):
            for j in range(len(fp[i])):
                fp[i][j] = fp[i][j] * (params['delta'] / params['utau'])
                Sp[i][j] = Sp[i][j] / ((params['rhom']*(params['utau'])**2)**2) 
        pl.plot_psd_pressure(fp, Sp, params['path'])
        # Velocity
        umag = objPSD.GetUmag()
        (fu, Su) = fc.CalculatePSD(umag, dt, params)
        for i in range(len(fu)):
            for j in range(len(fu[i])):
                fu[i][j] = fu[i][j] * (params['delta'] / params['utau'])
                Su[i][j] = Su[i][j] / (params['utau']**2)
        pl.plot_psd_velocity(fu, Su, params['path'])  
        # delete object  
        del objPSD

    if params['routine']['corr']:
        objCorr = DataSpatialCorrelation(params)
        up = objCorr.GetUprime()
        vp = objCorr.GetVprime()
        wp = objCorr.GetWprime()

        Rii = []
        Rii.append(fc.CalculateTwoPointsCorrelation(up[0], params))
        Rii.append(fc.CalculateTwoPointsCorrelation(vp[0], params))
        Rii.append(fc.CalculateTwoPointsCorrelation(wp[0], params))
        z = np.linspace(0, math.pi/2.0, num=params['npcorr']) / params['lzRef']
        pl.plot_tpcorr(z, Rii, params['path'])
        
        #Rii = []
        #Rii.append(fc.CalculateOnePointCorrelation(up[0], params))
        #Rii.append(fc.CalculateOnePointCorrelation(vp[0], params))
        #Rii.append(fc.CalculateOnePointCorrelation(wp[0], params))
        #t = np.linspace(0, params['time'], num=int(params['numSteps'] / 
        #                                        params['pfreq']) + 1) / params['tRef']
        #pl.plot_opcorr(t, Rii, params['path'])

        # delete object  
        del objCorr

    if params['routine']['profile']:
        # Data
        objProfE = DataDNSDataProfile(params)
        objProfEP = DataDNSDataProfilePressure(params)
        objRSE = DataDNSDataReyStress(params)
        objProf = DataTBLVProfile(params)

        yplus = []
        step = 0.0
        dy = 1.0/(params['npointsy'] - 1) 
        for i in range(params['npointsy']):
            yplus.append(step*(params['utau']/params['nu']))   
            step += dy

        #ypluspe : profile
        ypluspe = objProfE.GetYplus()

        # Uplus
        Uplus = objProf.GetUplus()
        Upluse = objProfE.GetUmean()
        headers = ['yplus', 'Uplus0', 'Uplus1', 'Uplus2']
        fc.WriteVectorsToCSV(params['path']+'/tblvprofile.csv', headers, yplus, Uplus)
        pl.plot_tblvprofile(yplus, Uplus, ypluspe, Upluse, params['path'])

        #yplure : reystresses
        yplusre = objRSE.GetYplus()

        #u_{rms}^{plus}
        urmsp = objProf.GetURMSP()
        urmspe = objRSE.GetURMS()
        headers = ['yplus', 'urms0', 'urms1', 'urms2']
        fc.WriteVectorsToCSV(params['path']+'/urms.csv', headers, yplus, urmsp)
        pl.plot_urmsp(yplus, urmsp, yplusre, urmspe, params['path'])
    
        #v_{rms}^{plus}    
        vrmsp = objProf.GetVRMSP()
        vrmspe = objRSE.GetVRMS()
        headers = ['yplus', 'vrms0', 'vrms1', 'vrms2']
        fc.WriteVectorsToCSV(params['path']+'/vrms.csv', headers, yplus, vrmsp)
        pl.plot_vrmsp(yplus, vrmsp, yplusre, vrmspe, params['path'])

        #w_{rms}^{plus}    
        wrmsp = objProf.GetWRMSP()
        wrmspe = objRSE.GetWRMS()
        headers = ['yplus', 'wrms0', 'wrms1', 'wrms2']
        fc.WriteVectorsToCSV(params['path']+'/wrms.csv', headers, yplus, wrmsp)
        pl.plot_wrmsp(yplus, wrmsp, yplusre, wrmspe, params['path'])

        # -uv / u_{tau}^{2}    
        uvp = objProf.GetUVP()
        uvpe = objRSE.GetUV()
        headers = ['yplus', 'uv0', 'uv1', 'uv2']
        fc.WriteVectorsToCSV(params['path']+'/uv.csv', headers, yplus, uvp)
        pl.plot_uvp(yplus, uvp, yplusre, uvpe, params['path'])

        # p_{rms}     
        prmsp = objProf.GetPRMSP()
        ppe = objProfEP.GetPRMSP()
        ypluse = objProfEP.GetYplus() # both
        headers = ['yplus', 'prms0', 'prms1', 'prms2']
        fc.WriteVectorsToCSV(params['path']+'/prms.csv', headers, yplus, prmsp)
        pl.plot_prmsp(yplus, prmsp, ypluse, ppe, params['path'])

        # up  
        upp = objProf.GetUPP()
        upe = objProfEP.GetUPP()
        headers = ['yplus', 'up0', 'up1', 'up2']
        fc.WriteVectorsToCSV(params['path']+'/up.csv', headers, yplus, upp)
        pl.plot_upp(yplus, upp, ypluse, upe, params['path'])

        # Analytical solution
        #uplusAnalyt = []
        #yplusAnalyt = []
        #y = -params['Ly']
        #while(y < 0):
        #    uplusAnalyt.append((1.36346011e+05*pow(1-abs(y), 14)-1.00157155e+06*pow(1-abs(y), 13)+3.30835371e+06*pow(1-abs(y), 12)-6.49079211e+06*pow(1-abs(y), 11)+8.41556333e+06*pow(1-abs(y), 10)-7.58969284e+06*pow(1-abs(y), 9)+4.87934144e+06*pow(1-abs(y), 8)-2.25289875e+06*pow(1-abs(y), 7)+7.41903706e+05*pow(1-abs(y), 6)-1.70090684e+05*pow(1-abs(y), 5)+2.56967948e+04*pow(1-abs(y), 4)-2.21227945e+03*pow(1-abs(y), 3)+4.33607541e+01*pow(1-abs(y), 2)+1.10361842e+01*pow(1-abs(y), 1)+3.28951994e-04*pow(1-abs(y), 0))/params['utau'])
        #    yplusAnalyt.append((y+1)*(params['utau']/params['nu'])) 
        #    y+=0.0001

        del objProf

else:
    # Experimental Data
    objProfE = DataDNSDataProfile(params)
    objRSE = DataDNSDataReyStress(params)
    objProfEP = DataDNSDataProfilePressure(params)

    #ypluspe : profile
    ypluspe = objProfE.GetYplus()

    # Uplus
    Ufplus = fc.LoadCSV(params['path']+'/tblvprofile.csv')
    Uplus = []
    Uplus.append(Ufplus['Uplus0'])
    Uplus.append(Ufplus['Uplus1'])
    Uplus.append(Ufplus['Uplus2'])
    Upluse = objProfE.GetUmean()
    pl.plot_tblvprofile(Ufplus['yplus'], Uplus, ypluspe, Upluse, params['path'])

    #yplure : reystresses
    yplusre = objRSE.GetYplus()

    #u_[rms]^{plus}
    urmspe = objRSE.GetURMS()
    ufrms = fc.LoadCSV(params['path']+'/urms.csv')
    urms = []
    urms.append(ufrms['urms0'])
    urms.append(ufrms['urms1'])
    urms.append(ufrms['urms2'])
    pl.plot_urmsp(ufrms['yplus'], urms, yplusre, urmspe, params['path'])

    #v_[rms]^{plus}    
    vrmspe = objRSE.GetVRMS()
    vfrms = fc.LoadCSV(params['path']+'/vrms.csv')
    vrms = []
    vrms.append(vfrms['vrms0'])
    vrms.append(vfrms['vrms1'])
    vrms.append(vfrms['vrms2'])
    pl.plot_vrmsp(vfrms['yplus'], vrms, yplusre, vrmspe, params['path'])

    #w_[rms]^{plus}    
    wrmspe = objRSE.GetWRMS()
    wfrms = fc.LoadCSV(params['path']+'/wrms.csv')
    wrms = []
    wrms.append(wfrms['wrms0'])
    wrms.append(wfrms['wrms1'])
    wrms.append(wfrms['wrms2'])
    pl.plot_wrmsp(wfrms['yplus'], wrms, yplusre, wrmspe, params['path'])

    # -uv / u_{tau}    
    uvpe = objRSE.GetUV()
    uvf = fc.LoadCSV(params['path']+'/uv.csv')
    uv = []
    uv.append(uvf['uv0'])
    uv.append(uvf['uv1'])
    uv.append(uvf['uv2'])
    pl.plot_uvp(uvf['yplus'], uv, yplusre, uvpe, params['path'])
        
    # p_{rms}^{plus}    
    ppe = objProfEP.GetPRMSP()
    ppf = fc.LoadCSV(params['path']+'/prms.csv')
    pp = []
    pp.append(ppf['prms0'])
    pp.append(ppf['prms1'])
    pp.append(ppf['prms2'])
    pl.plot_prmsp(ppf['yplus'], pp, yplusre, ppe, params['path'])

    # up     
    upe = objProfEP.GetUPP()
    upf = fc.LoadCSV(params['path']+'/up.csv')
    up = []
    up.append(upf['up0'])
    up.append(upf['up1'])
    up.append(upf['up2'])
    pl.plot_upp(upf['yplus'], up, yplusre, upe, params['path'])

    
