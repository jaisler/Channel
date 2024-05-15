import matplotlib.pyplot as plt
from matplotlib import rc, cm
import numpy as np
import math
rc('text', usetex=True)

def plot_psd_pressure(f,S,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})

    # -5/3 power law
    #x0=1e4
    #x1=1e5
    #x = np.linspace(x0, x1)
    #y = -(5/3)*10*np.log10(x) +(5/3)*10*np.log10(x0) 
    #p0, = plt.semilogx(x, y, '-', color='k')
    
    p0, = plt.loglog(f[0], S[0], '-', color='r', linewidth=3)	
    p1, = plt.loglog(f[1], S[1], '-', color='b', linewidth=3)		
    p2, = plt.loglog(f[2], S[2], '-', color='g', linewidth=3)
    
    
    plt.legend([p0,p1,p2],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.123, right=0.96, bottom=0.135, top=0.96)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.xlim((1e-1,500))
    plt.ylim((1e-6, 10))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$f \delta / u_{\tau}$',fontsize = 18)
    plt.ylabel(r'$PSD(p) / (\rho u_{\tau}^{2})^{2}$',fontsize = 18)
    fig.savefig(path + '/psdp.pdf', format='PDF')
    fig.savefig(path + '/psdp.png', format='png')
    plt.show()

def plot_psd_velocity(f,S,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})

    p0, = plt.loglog(f[0], S[0], '-', color='r', linewidth=3)	
    p1, = plt.loglog(f[1], S[1], '-', color='b', linewidth=3)		
    p2, = plt.loglog(f[2], S[2], '-', color='g', linewidth=3)	
    
    plt.legend([p0,p1,p2],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((1e-1,500))
    plt.ylim((1e-6, 5))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$f \delta / u_{\tau}$',fontsize = 18)
    plt.ylabel(r'$PSD(u) / u_{\tau}^{2}$',fontsize = 18)
    fig.savefig(path + '/psdu.pdf', format='PDF')
    fig.savefig(path + '/psdu.png', format='png')
    plt.show()

def plot_tblvprofile(y,u,ya,ua,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ya, ua, '-', color='k', linewidth=2)
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    #plt.ylim((1e-11, 0.2))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$U^{+}$',fontsize = 18)
    fig.savefig(path + '/prof.pdf', format='PDF')
    fig.savefig(path + '/prof.png', format='png')
    plt.show()

def plot_urmsp(y,u,ye,ue,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    #plt.ylim((-0.1, 3.0))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$u_{rms}^{+}$',fontsize = 18)
    fig.savefig(path + '/prof_urmsp.pdf', format='PDF')
    fig.savefig(path + '/prof_urmsp.png', format='png')
    plt.show()

def plot_vrmsp(y,u,ye,ue,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')	
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    #plt.ylim((-0.1, 1.0))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$v_{rms}^{+}$',fontsize = 18)
    fig.savefig(path + '/prof_vrmsp.pdf', format='PDF')
    fig.savefig(path + '/prof_vrmsp.png', format='png')
    plt.show()

def plot_wrmsp(y,u,ye,ue,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')	
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    #plt.ylim((-0.1, 1.5))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$w_{rms}^{+}$',fontsize = 18)
    fig.savefig(path + '/prof_wrmsp.pdf', format='PDF')
    fig.savefig(path + '/prof_wrmsp.png', format='png')
    plt.show()

def plot_uvp(y,u,ye,ue,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')	
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    #plt.ylim((-0.1, 1.0))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$-\langle uv \rangle / u_{\tau}^{2}$',fontsize = 18)
    fig.savefig(path + '/prof_uvp.pdf', format='PDF')
    fig.savefig(path + '/prof_uvp.png', format='png')
    plt.show()

def plot_prmsp(y,u,ye,ue,path): 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')	
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe],
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    plt.ylim((-0.1, 4.0))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$p_{rms}^{+}$',fontsize = 18)
    fig.savefig(path + '/prof_prmsp.pdf', format='PDF')
    fig.savefig(path + '/prof_prmsp.png', format='png')
    plt.show()

def plot_upp(y,u,ye,ue,path): 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})
    plt.xscale('log')

    pe, = plt.plot(ye, ue, '-', color='k', linewidth='2')	
    p0, = plt.plot(y, u[0], 'o', color='r', mfc='none', markersize='7')	
    p1, = plt.plot(y, u[1], 's', color='b', mfc='none', markersize='7')
    p2, = plt.plot(y, u[2], '^', color='g', mfc='none', markersize='7')	
    
    plt.legend([p0,p1,p2,pe], #,pe
    [r'$x/\delta=5$',
     r'$x/\delta=10$',
     r'$x/\delta=55$',
     r'DNS data'],
     loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((6e-1,150))
    plt.ylim((-0.1, 0.8))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{+}$',fontsize = 18)
    plt.ylabel(r'$\langle up \rangle / \rho u_{\tau}^{3}$',fontsize = 18)
    fig.savefig(path + '/prof_upp.pdf', format='PDF')
    fig.savefig(path + '/prof_upp.png', format='png')
    plt.show()

def plot_tpcorr(z, R, path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})

    p0, = plt.plot(z, R[0], 'o-', color='b', markersize='6')	
    p1, = plt.plot(z, R[1], 's-', color='r', markersize='6')	
    p2, = plt.plot(z, R[2], '^-', color='g', markersize='6')	

    plt.legend([p0, p1, p2],
    [r'$R_{u^{\prime}u^{\prime}}$',
     r'$R_{v^{\prime}v^{\prime}}$',
     r'$R_{w^{\prime}w^{\prime}}$'],
    # r'DNS data'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((0,6.3))
    plt.ylim((-0.1, 1.1))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$z / l_{z}^{ref}$',fontsize = 18)
    plt.ylabel(r'$R_{i^{\prime}i^{\prime}}$',fontsize = 18)
    fig.savefig(path + '/corr.pdf', format='PDF')
    fig.savefig(path + '/corr.png', format='png')
    plt.show()

def plot_opcorr(z, R, path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})

    p0, = plt.plot(z, R[0], 'o-', color='b', markersize='6')	
    #p1, = plt.plot(z, R[1], 's-', color='r', markersize='6')	
    #p2, = plt.plot(z, R[2], '^-', color='g', markersize='6')	

    #plt.legend([p0, p1, p2],
    #[r'$R_{u^{\prime}u^{\prime}}$',
    # r'$R_{v^{\prime}v^{\prime}}$',
    # r'$R_{w^{\prime}w^{\prime}}$'],
    # r'DNS data'],
    #loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((0,4))
    #plt.ylim((-0.1, 1.1))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$t / t^{ref}$',fontsize = 18)
    plt.ylabel(r'$R_{i^{\prime}i^{\prime}}$',fontsize = 18)
    fig.savefig(path + '/corr.pdf', format='PDF')
    fig.savefig(path + '/corr.png', format='png')
    plt.show()