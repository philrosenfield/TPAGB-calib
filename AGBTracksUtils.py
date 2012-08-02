import trilegal_diagnostics
import GoogleSitesTable
import galaxy_tests
import os,sys
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=False
rcParams['axes.linewidth'] = 1
rcParams['ytick.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['axes.edgecolor'] = 'black'
#rc('text',usetex=True)
from matplotlib.ticker import *
nullfmt = NullFormatter() # no labels

def discrete_colors(Ncolors,colormap='gist_rainbow'):
    colors = []
    cmap = cm.get_cmap(colormap)
    for i in range(Ncolors):
        colors.append( cmap(1.*i/Ncolors) ) # color will now be an RGBA tuple
    return colors

def load_input(filename):
    '''
    reads an input file into a dictionary. 
    file must have key first then value(s)
    Will make 'True' into a boolean True
    Will understand if a value is a float, string, or list, etc.
    Ignores all lines that start with #, but not with # on the same line as key, value.
    '''
    try:
        literal_eval
    except NameError:
        from ast import literal_eval
    
    lines = open(filename).readlines()
    d = {}
    for line in lines:
        if line.startswith('#'): continue
        if len(line.strip()) == 0: continue
        key,val = line.strip().partition(' ')[0::2]
        d[key] = is_numeric(val.replace(' ',''))
    # do we have a list?
    for key in d.keys():
        # float
        if type(d[key]) == float: continue
        # list:
        temp = d[key].split(',')
        if len(temp) > 1:
            try:
                d[key] = map(float,temp)    
            except:
                d[key] = temp
        # dict:
        elif len(d[key].split(':'))>1:
            temp1 = d[key].split(':')
            d[key]={is_numeric(temp1[0]):is_numeric(temp1[1])}
        else:        
            val = temp[0]
            # boolean
            if val.upper().startswith('TR') or val.upper().startswith('FA'): val = literal_eval(val)
            # string
            d[key] = val
    return d

class AGBTracks(object):
    '''
    A fast way of loading files.
    
    examples:
    track = some agb track file name (string)
    AGB = get_numeric_data(track)
    
    all the data: 
    AGB.data_array
    
    the file name:
    AGB.name
    
    a dictionary with 'key': column_number
    AGB.key_dict
    
    get the row of the AGB track by calling the number in the first column (ex: 7):
    AGB.get_row_bynum(7)
    
    get a specified column (use AGB.key_dict to know which one, or do head -1) (ex: logL)
    AGB.get_col('L_*')
    '''
    def __init__(self, data_array, col_keys,name):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys,range(len(col_keys))))
        self.name = name
        self.mass = float(os.path.split(name)[1].split('_')[1])
        self.metallicity = float(os.path.split(name)[1].split('_')[2].replace('Z',''))
        
    def get_row(self,i):
        return self.data_array[i,:]
    
    def get_row_bynum(self,i):
        row = np.nonzero(self.data_array[:,self.key_dict['step']] == i)[0]
        return self.data_array[row,:]
    
    def get_col(self,key):
        return self.data_array[key]

def get_numeric_data(filename):
    '''
    made to read all of Paola's tracks. It takes away her "lg" meaning log.
    Returns an AGBTracks object. If there is a problem reading the data, all data
    are passed as zeros.
    '''
    f=open(filename,'r')
    col_keys = f.readline().replace('#','').replace('lg','').replace('*','star').strip().split()
    f.close()
    try:
        data = np.genfromtxt(filename,missing_values='************',names=col_keys)
    except ValueError:
        print 'problem with',filename
        data = np.zeros(len(col_keys))
    return AGBTracks(data,col_keys,filename)
        
def make_iso_file(track,Qs,slopes,isofile):
    '''
    this only writes the quiescent lines and the first line.
    format of this file is:
    t_min,          age in yr
    logl_min,       logL
    logte_min,      logTe
    mass_min,       actual mass along track
    mcore_min,      core mass
    co_min,         C/O ratio
    per_min,        period in days
    ip_min,         1=first overtone, 0=fundamental mode
    mlr_min,        - mass loss rate in Msun/yr
    logtem_min,     keep equal to logTe
    x_min,          X
    y_min,          Y
    xcno_min        X_C+X_O+X_N
    slope           dTe/dL
    '''
    fmt = '%.4e %.4f %.4f %.5f %.5f %.4f %.4e %i %.4e %.4f %.6e %.6e %.6e %.4f\n'

    # cull agb track to quiescent, write out.
    #rows = [q-1 for q in Qs[:-1]] # cutting the final point for the file.
    rows = [q-1 for q in Qs] # cutting the final point for the file.
    rows[0] +=2
    
    #f=open(track.name,'r').readlines()
    #col_keys = f[0].replace('#','').replace('lg','').strip().split()
    col_keys = np.array(track.key_dict.keys())[np.argsort(track.key_dict.values())]
    #check on P0 == per_min
    #outcols = ['ageyr','L_star','T_star','M_star','M_c','C/O','P0','P1','Pmod','dMdt','T_star','H','Y']
    #outcols.extend([key for key in col_keys if (key.startswith('C1') or key.startswith('N1') or key.startswith('O1'))])
    cno = [key for key in col_keys if (key.startswith('C1') or key.startswith('N1') or key.startswith('O1'))]
    #col_nums = [int(list(col_keys).index(key)) for key in outcols]
    isofile.write('%.4f %i # %s\n'%(track.mass,len(rows),os.path.split(track.name)[1]))
    for r in rows:
        try:
            #data_row = map(float,f[r].split())
            data_row = track.data_array[r]
        except ValueError:
            print r
            print ''
            print rows
            print ''
            print f[r-1]
            print track.name
            return
        CNO = np.sum([data_row[c] for c in cno])
        mdot = 10**(data_row['dMdt'])
        if data_row['Pmod'] == 0:
            period = data_row['P0']
        else:
            period = data_row['P1']
        if r == rows[-1]: # adding nonsense slope for the final row.
            slope = 999999
        else:
            slope = 1./slopes[list(rows).index(r)]
        try:
            isofile.write(fmt%(data_row['ageyr'],data_row['L_star'],data_row['T_star'],
                               data_row['M_star'],data_row['M_c'],data_row['CO'],period,
                               data_row['Pmod'],mdot,data_row['T_star'],data_row['H'],
                               data_row['Y'],CNO,slope))
                               

            #isofile.write('%.3f %.4f %.4f %.5f %.5f %.4f %.4e %i %.4e %.4f %.6e %.6e %.6e %.4f\n'%(data_row[col_nums[0]],data_row[col_nums[1]],data_row[col_nums[2]],data_row[col_nums[3]],data_row[col_nums[4]],data_row[col_nums[5]],period,data_row[col_nums[8]],mdot,data_row[col_nums[2]],data_row[col_nums[11]],data_row[col_nums[12]],CNO,1./slopes[list(rows).index(r)]))
        except IndexError:
            print list(rows).index(r)
            print len(rows),len(slopes)
            print 1./slopes[list(rows).index(r)]
    return

def write_cmd_input_file(**kwargs):
    '''
    make a TRILEGAL cmd_input file based on default. 
    
    Send each parameter that is different than default by:
    kwargs = { 'kind_tpagb': 4, 'file_tpagb': 'isotrack/tracce_CAF09_S0.dat'}
    cmd_input_file = write_cmd_input_file(**kwargs)   
    
    To make the default file:
    cmd_input_file = write_cmd_input_file()
    
    if you don't specify cmd_input_file, output goes to cmd_input_TEMP.dat
    '''
    kind_tracks = kwargs.get('kind_tracks',2)
    file_isotrack = kwargs.get('file_isotrack','isotrack/parsec/CAF09.dat')
    file_lowzams = kwargs.get('file_lowzams','isotrack/bassazams_fasulla.dat')
    kind_tpagb = kwargs.get('kind_tpagb',4)
    file_tpagb = kwargs.get('file_tpagb','isotrack/isotrack_agb/tracce_CAF09_AFEP02_I1_S1.dat')
    kind_postagb = kwargs.get('kind_postagb',0)
    file_postagb = kwargs.get('file_postagb','isotrack/final/pne_wd_test.dat')
    mass_loss = kwargs.get('mass_loss')
    if mass_loss:
        kind_rgbmloss = 1
        law_mass_loss, = mass_loss.keys()
        efficiency_mass_loss, = mass_loss.values()
    # these are for using cmd2.2:
    kind_mag = kwargs.get('kind_mag',None)
    photsys = kwargs.get('photsys','wfpc2')
    file_mag = 'tab_mag_odfnew/tab_mag_'+photsys+'.dat'
    kind_imf = kwargs.get('kind_imf',None)
    file_imf = kwargs.get('file_imf','tab_imf/imf_chabrier_lognormal.dat')
    
    # if not using cmd2.2:
    if kind_imf == None:
        kind_imfr = kwargs.get('kind_imfr',0)
        file_imfr = kwargs.get('file_imfr','tab_ifmr/weidemann.dat')
    
    track_comments = '# kind_tracks, file_isotrack, file_lowzams'
    tpagb_comments = '# kind_tpagb, file_tpagb'
    pagb_comments = '# kind_postagb, file_postagb DA VERIFICARE file_postagb'
    mag_comments = '# kind_mag, file_mag'
    imf_comments = '# kind_imf, file_imf'
    imfr_comments = '# ifmr_kind, file with ifmr'
    mass_loss_comments = '# RGB mass loss: kind_rgbmloss, law, and its efficiency'
    footer = '################################explanation######################\nkind_tracks: 1= normal file \n file_isotrack: tracks for low+int mass\nfile_lowzams: tracks for low-ZAMS \n kind_tpagb: 0= none\n        1= Girardi et al., synthetic on the flight, no dredge up \n        2= Marigo & Girardi 2001, from file, includes mcore and C/O\n        3= Marigo & Girardi 2007, from file, includes period, mode and mloss  \n        4= Marigo et al. 2011, from file, includes slope  \n file_tpagb: tracks for TP-AGB\n\nkind_postagb: 0= none\n          1= from file\nfile_postagb: PN+WD tracks\n\nkind_ifmr: 0= default\n           1= from file\n'
    cmd_input_file = kwargs.get('cmd_input_file','cmd_input_TEMP.dat')
    fh = open(cmd_input_file,'w')
    formatter = '%i %s %s \n'
    fh.write('%i %s %s %s \n' % (kind_tracks,file_isotrack,file_lowzams,track_comments))
    fh.write(formatter % (kind_tpagb,file_tpagb,tpagb_comments))
    fh.write(formatter % (kind_postagb,file_postagb,pagb_comments))
    if kind_mag != None: fh.write(formatter % (kind_mag,file_mag,mag_comments))
    if kind_imf == None: 
        fh.write(formatter % (kind_imfr,file_imfr,imfr_comments))
    else:
        fh.write(formatter % (kind_imf,file_imf,imf_comments))
    if mass_loss: 
        fh.write('%i %s %.3f \n'%(kind_rgbmloss,law_mass_loss,efficiency_mass_loss))
    fh.write(footer)
    fh.close()
    return cmd_input_file

def make_met_file(tracce,Zs,Ys,isofiles):
    t = open(tracce,'w')
    t.write('%i\n'%len(isofiles))
    [t.write('%.4f\t%.3f\t%s\n'%(Zs[i],Ys[i],isofiles[i])) for i in np.argsort(Zs)]
    t.close()
    return

def get_TP_inds(ntp):
    un =  np.unique(ntp)
    iTPs = [list(ntp).index(u) for u in un] # this is the first step in each TP.             
    TPs = [np.arange(iTPs[i],iTPs[i+1]) for i in range(len(iTPs)-1)] # The indices of each TP.
    TPs.append(np.arange(iTPs[i+1],len(ntp))) # don't forget the last one.
    return TPs

def get_unique_inds(ntp):
    un =  np.unique(ntp)
    iTPs = [list(ntp).index(u) for u in un] # this is the first step in each TP.             
    TPs = [np.arange(iTPs[i],iTPs[i+1]) for i in range(len(iTPs)-1)] # The indices of each TP.
    TPs.append(np.arange(iTPs[i+1],len(ntp))) # don't forget the last one.
    return TPs,iTPs

def add_points_to_q_track(track,qs):
    '''
    when to add an extra point for low masses
    if logt[qs+1] is hotter than logt[qs]
    and there is a point inbetween logt[qs] and logt[qs+1] that is cooler
    than logt[qs]
    add the coolest point.
    '''
    addpt=[]
    logt = AGBTracks.get_col(track,'T_star')
    step = AGBTracks.get_col(track,'step')
    status =AGBTracks.get_col(track,'status')
    Tqs = logt[qs]
    # need to use some unique array, not log t, since log t could repeat,
    # index would find the first one, not necessarily the correct one.
    Sqs = step[qs]-1. # makes it the same as qs.
    # takes the difference in logt(qs) to see if we get hotter or colder.
    # finds where the logt goes from getting colder to hotter...
    ht = np.where(np.diff(np.sign(np.diff(Tqs))))[0]+1
    ht = np.append(ht,ht+1) # between the first and second
    Sqs_ht = Sqs[ht]
    # the indices between each hot point.
    t_mids = [map(int,step[int(Sqs_ht[i]):int(Sqs_ht[i+1])]) for i in range(len(Sqs_ht)-1)]
    Sqs_ht = Sqs_ht[:-1]
    for i in range(len(Sqs_ht)):
        hot_inds = np.nonzero(logt[int(Sqs_ht[i])] > logt[t_mids[i]])[0]
        if len(hot_inds) > 0:
            # index of the min T of the hot index from above.
            addpt.append(list(logt).index(np.min(logt[[t_mids[i][hot_ind] for hot_ind in hot_inds]])))
            
    if len(addpt)>0: addpt = [a for a in addpt if status[a]==7.]
    Qs = np.sort(np.concatenate((addpt,qs)))
    return Qs,addpt

def find_dldt(track,TPs,addpt):
    '''
    Finds dL/dt of track object
    '''
    status = AGBTracks.get_col(track,'status')
    #dl/dT seems somewhat linear for 0.2<phi<0.4 ...
    phi = AGBTracks.get_col(track,'PHI_TP')
    phi[0]= -1. # The first line in the agb track is 1. This isn't the quiessent...
    lin_rise = np.nonzero((status==7) & (phi<0.4) & (phi>0.2))[0]
    rising = [list(set(TP) & set(lin_rise)) for TP in TPs] 
    logl = AGBTracks.get_col(track,'L_star')
    logt = AGBTracks.get_col(track,'T_star')
    order = 1
    fits = [np.polyfit(logt[r],logl[r],order) for r in rising if len(r)>0]
    slopes = [fits[i][0] for i in range(len(fits))]
    if len(addpt)>0:
        addrise = [TPs.index(TP) for TP in TPs if len(set(addpt) & set(TP)) > 0 ]
        addslope = slopes[addrise[0]]
        Slopes = []
        for s in slopes:
            if s == addslope:
                Slopes.append(addslope)
            Slopes.append(s)
    else:
        Slopes = slopes[:]
    return rising,Slopes,fits

def set_up_three_panel_plot():
    fig = plt.figure(figsize=(8,8))
    
    d=0.02
    low_b = 0.1 
    
    left, width= 0.1,0.8
    dh = (1.-3*d-low_b-.05)/3.
    
    low_t = low_b + dh
    mid_b = low_t + d
    mid_t = mid_b + dh
    top_b = mid_t + d
    top_t = top_b + dh
    
    bottoms = [low_b, mid_b, top_b]
    
    axs = [plt.axes([left,b,width,dh]) for b in bottoms]
    return axs

def two_panel_plot_vert():
    fig = plt.figure(2, figsize=(8,8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.41
    dh = 0.03
    
    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom+height+dh), width, height]
    
    ax1 = plt.axes(axis1)
    ax2 = plt.axes(axis2)
    ax2.xaxis.set_major_formatter(nullfmt)
    
    return ax1,ax2

def diag_plots(track,logl,logt,age,slopes,Qs,addpt,rising,fits,**kwargs):
    agb_mix = kwargs.get('agb_mix')
    set_name = kwargs.get('set_name')
    ext = '.png'
    logt_lim = (3.75,3.35)
    logl_lim = (2.4,4.8)
    lage_lim = (1.,10**7)
    co_lim = (0,5)
    #x = np.linspace(min(logt),max(logt))
    
    # HRD
    fig = plt.figure()
    ax = plt.axes()
    plotpath = os.path.join(kwargs.get('diagnostic_dir'),'HRD/')
    ensure_dir(plotpath)
    ax.annotate(r'$%s$'%agb_mix.replace('_','\ '),xy=(3.43,2.8))
    ax.annotate(r'$%s$'%set_name,xy=(3.43,2.7))
    ax.annotate(r'$M=%.2f$'%track.mass,xy=(3.43,2.6))
    #[plt.plot(x,fits[i][0]*x+fits[i][1],'--',color='blue',lw=0.5) for i in range(len(fits))]
    ax.plot(logt,logl,color='black')
    #[plt.plot(logt[r],logl[r],color='red',lw=2) for r in rising if len(r)>0]
    ax.plot(logt[map(int,Qs)],logl[map(int,Qs)],color='green',lw=2)
    [ax.plot(logt[q],logl[q],'o',color='green') for q in Qs]
    [ax.plot(logt[add],logl[add],'o',color='purple') for add in addpt if len(addpt)>0]
    ax.set_xlim(logt_lim)
    ax.set_ylim(logl_lim)
    ax.set_xlabel(r'$\log\ Te$')
    ax.set_ylabel(r'$\log\ L_{\odot}$')
    fig_name = os.path.join(plotpath,'_'.join(('diag',os.path.split(track.name)[1].replace('.dat',''))))
    plt.savefig('%s%s'%(fig_name,ext))
    #plt.axis([max(logt)+0.01,min(logt)-0.01,min(logl)-0.01,max(logl)+0.01])
    plt.close()
    
    ax1,ax2,ax3 = set_up_three_panel_plot()
    # AGE v TE
    ax1.plot(age,logt,color='black')
    ax1.plot(age[Qs],logt[Qs],'o',color='green')
    ax1.yaxis.set_major_locator(MultipleLocator(.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(.05))
    # AGE v L
    ax2.plot(age,logl,color='black')
    ax2.plot(age[Qs],logl[Qs],'o',color='green')
    ax2.yaxis.set_major_locator(MultipleLocator(.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(.1))
    # AGE vs CO
    CO = track.get_col('CO')
    ax3.plot(age,CO,color='black')
    ax3.plot(age[Qs],CO[Qs],'o',color='green')
    majorFormatter = ScalarFormatter()
    majorFormatter.set_powerlimits((-3, 4))
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.annotate(r'\rm{%s}'%agb_mix.replace('_','\ '),xy=(.06,.87),textcoords =('axes fraction'))
    ax3.annotate(r'\rm{%s}'%set_name.replace('_','\ '),xy=(.06,.77),textcoords =('axes fraction'))
    ax3.annotate(r'$M=%.2f$'%track.mass,xy=(.06,.67), textcoords =('axes fraction') )
    ax3.set_ylabel(r'$C/O$')
    ax2.set_ylabel(r'$\log\ L_{\odot}$')
    ax1.set_xlabel(r'$\rm{Age (yr)}$')
    ax1.set_ylabel(r'$\log\ Te$')
    #ax1.set_ylim(logt_lim)
    #ax2.set_ylim(logl_lim)
    #ax3.set_ylim(co_lim)
    
    plotpath = os.path.join(kwargs.get('diagnostic_dir'),'age_v/')
    ensure_dir(plotpath)
    fig_name = os.path.join(plotpath,'_'.join(('diag',os.path.split(track.name)[1].replace('.dat',''))))

    majorLocator = MultipleLocator()
    for ax in ax1,ax2,ax3:
        #ax.set_xscale('log')
        #ax.yaxis.set_major_locator(majorLocator)  
        ax.grid(linestyle=':',color='grey')
    plt.savefig('%s%s%s'%(fig_name,'_age_v',ext))
    plt.close()
    
    

    return

def make_readme(agb_iso_track,agb_mix):
    readmef = agb_iso_track+'readme.txt'
    if os.path.isfile(readmef): 
        readme = open(agb_iso_track+'readme.txt','a')
        readme.write('\n --- %s --- \n'%time.strftime("%Y-%m-%d %H:%M:%S"))
    else:
        readme = open(agb_iso_track+'readme.txt','w')
    readme.write('This directory has agb tracks parsed from AGBTracksUtils. Sources are \n %s \n and within\n'%agbtrack_dir)
    readme.close()

def plot_ifmr(imfrfile):
    mi,mf,z = np.loadtxt(imfrfile,unpack=True)
    zinds,unz = get_unique_inds(z)
    cols = discrete_colors(len(zinds))
    [plt.plot(mi[zinds[i]],mf[zinds[i]],color=cols[i],label=str(z[unz[i]])) for i in range(len(unz))]
    [plt.plot(mi[zinds[i]],mf[zinds[i]],'o',ms=4,color=cols[i],mec='white',alpha=0.5) for i in range(len(unz))]
    plt.legend(loc=2,frameon=False)
    plt.xlabel(r'$M_i/M_{\odot}$',fontsize=20)
    plt.ylabel(r'$M_f/M_{\odot}$',fontsize=20)
    plt.savefig(imfrfile.replace('dat','png'))
    plt.close()
    return

def two_panel_plot_vert(fign=2):
    fig = plt.figure(fign, figsize=(8,8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.41
    dh = 0.03
    
    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom+height+dh), width, height]
    
    ax1 = plt.axes(axis1)
    ax2 = plt.axes(axis2)
    ax2.xaxis.set_major_formatter(nullfmt)
    return ax1,ax2

def plot_cluster_test(lifetimesfile,**kwargs):
    agbtrack_dir = kwargs.get('agbtrack_dir')
    cluster_data = os.path.join(agbtrack_dir,'cluster_data.dat')
    cmass, ctauc, ecp, ecm, ctaum, emp, emm = np.loadtxt(cluster_data,unpack=True)
    z, mass, tauc, taum = np.loadtxt(lifetimesfile,unpack=True)
    zinds,unz = get_unique_inds(z)
    cols = discrete_colors(len(zinds))
        
    ax1,ax2 = two_panel_plot_vert(fign=3)
    ax2.set_ylabel(r'$\tau_M\ (\rm{Myr})$',fontsize=20)
    ax1.set_ylabel(r'$\tau_C\ (\rm{Myr})$',fontsize=20)
    ax1.set_xlabel(r'$M_{TO}\ (\rm{M}_\odot)$',fontsize=20)
    ax1.errorbar(cmass, ctauc, yerr=[ecm,ecp], fmt='o',color='black',label='LMC')
    ax2.errorbar(cmass, ctaum, yerr=[emm,emp], fmt='o',color='black')

    [ax1.plot(mass[zinds[i]],tauc[zinds[i]],color=cols[i],label=str(z[unz[i]])) \
     for i in range(len(unz))]
    ax1.legend(loc=0,frameon=False,numpoints=1)
    [ax2.plot(mass[zinds[i]],taum[zinds[i]],color=cols[i]) for i in range(len(unz))]
    
    plt.savefig(lifetimesfile.replace('dat','png'))
    print 'wrote',lifetimesfile.replace('dat','png')
    plt.close()
    return
    
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.isdir(d):
        os.makedirs(d)
        print 'made dirs:',d

def get_afile(src,search_string):
    '''
    simple search, returns a list
    '''
    import glob
    try:
        files = glob.glob1(src,search_string)
    except IndexError:
        print 'Can''t find',search_string,'in',src
        sys.exit(2)
    return [os.path.join(src,f) for f in files]
    
def is_numeric(lit):
    """
    Return value of numeric literal string or ValueError exception
    From http://rosettacode.org/wiki/Determine_if_a_string_is_numeric#Python
    """
    # Empty String
    if len(lit) <= 0:
        return lit    
    # Handle '0'
    if lit == '0': return 0
    # Hex/Binary
    if len(lit) > 1: # sometimes just '-' means no data...
        litneg = lit[1:] if lit[0] == '-' else lit
        if litneg[0] == '0':
            if litneg[1] in 'xX':
                return int(lit,16)
            elif litneg[1] in 'bB':
                return int(lit,2)
            else:           
                try:
                    return int(lit,8)
                except ValueError:
                    pass
    # Int/Float/Complex
    try:
        return int(lit)
    except ValueError:
        pass
    try:
        return float(lit)
    except ValueError:
        pass
    try:
        return complex(lit)
    except ValueError:
        pass
    return lit

def make_local_copy(file,dest=os.environ['ISOTRACK'][:-1]+'_agb/'):
    if dest!=None:
        os.system('cp %s %s'%(file,dest))
    return
    
def do_everything(**kwargs):
    '''
    HOW TO USE THIS:
    This script formats Paola's tracks and creates the files needed to use them with
    TRILEGAL.
    
    You can either pass a dictionary with the necessary locations of the files to be made
    and the files made or simply make an input file:
    
    kwargs = load_input(filename)
    do_everything(**kwargs)
    
    the input file should be formatted in this way (comments are optional):
    # Paola's tracks are here:
    agbtrack_dir    /Users/phil/research/TP-AGBcalib/AGBTracks    
    # This agb mix name
    agb_mix         agb_caf09
    # This set name
    set_name        S0
    # Where cmd_input files will go for trilegal
    trilegal_dir    output_for_trilegal/
    # Paola's formatted tracks are here:
    isotrack_dir    isotrack/parsec/AGB_TRACKS/
    # File to link from cmd_input to Paola's formatted tracks:
    tracce_dir      isotrack/parsec/
    # make initial and final mass relation? (Need more than one metallicity)
    make_imfr       False
    # Diagnostic plots will mimic directory structure of abgtrack_dir
    diagnostic_dir  DIAGPLOTS
    
    WHAT'S GOING ON:    
    Here is what it takes to go from Paola's agb tracks to trilegal:
    
    1. Paola's tracks
    2. Paola's tracks formatted for trilegal
    3. trilegal files that link to Paola's formatted tracks
    
    1. Paola's tracks should have naming scheme [mix]/[set]/[Trash]_[metallicity]_[Y]/[track_identifier]
    ex agb_caf09_z0.008/S1/agb_*Z*.dat
    ex CAF09/S_SCS/*Z*/agb_*Z*.dat
    These can go anywhere, and do not need to be included in trilegal directories.
    
    2. Paola's tracks are formatted for trilegal with no header file.
    These go in trilegal_1.3/isotrack/AGB_TRACKS/Z[metallicity]_[mix]_[set].dat
    They only include the quiescent phases and also a column with dL/dT
        
    3. trilegal needs two files to link to Paola's formatted tracks
       a. track file goes here: trilegal_1.3/isotrack/tracce_[mix]_[set].dat
       b. cmd_input_file that links to the track file
          goes here trilegal_1.3/cmd_input_[mix]_[set].dat
    '''
    home = os.getcwd()
    agbtrack_dir = kwargs.get('agbtrack_dir')
    agb_mix = kwargs.get('agb_mix')
    set_name = kwargs.get('set_name')
    trilegal_dir = kwargs.get('trilegal_dir')
    isotrack_dir = kwargs.get('isotrack_dir')
    tracce_dir  = kwargs.get('tracce_dir')
    diagnostic_dir0 = kwargs.get('diagnostic_dir')
    make_imfr = kwargs.get('make_imfr')
    make_copy = kwargs.get('make_copy')
    mass_loss = kwargs.get('mass_loss')
    # do dirs exist?
    ensure_dir(isotrack_dir)
    ensure_dir(os.path.join(diagnostic_dir0,agb_mix,set_name+'/'))
    # name convention: [mix]_[set].dat
    name_conv = '_'.join((agb_mix,set_name))+'.dat'
    
    if make_imfr == True:
        ifmr = os.path.join(diagnostic_dir0,agb_mix,set_name,'_'.join(('ifmr',name_conv)))
        mfile = open(ifmr,'w')
        mfile.write('# M_i M_f Z\n')
  
    # name convention: [mix]_[set].dat
    name_conv = '_'.join((agb_mix,set_name))+'.dat'
    
    # cmd_input_file that links to the track file
    cmd_input = '_'.join(('cmd','input',name_conv))
    cmd_input_file = os.path.join(home,trilegal_dir,cmd_input)
    
    # track file to link from cmd_input to paola's formatted tracks
    tracce_file = os.path.join(home,tracce_dir,'_'.join(('tracce',name_conv)))
    tracce_file_rel = os.path.join('isotrack',tracce_dir,'_'.join(('tracce',name_conv)))
    # to find AGB tracks (this could be a kwarg... 
    # I have Z because agb_imfr.dat would crash.
    track_identifier = 'agb_*Z*.dat'
    
    isofiles,Zs,Ys = [],[],[]
    
    # moving to the the directory with metallicities.
    working_dir = os.path.join(agbtrack_dir,agb_mix,set_name)
    os.chdir(working_dir)
    metal_dirs = [m for m in os.listdir(working_dir) if os.path.isdir(m)]
    metals = np.argsort([float(m.split('_')[1].replace('Z','')) for m in metal_dirs])
    print 'found %i metallicities'%len(metal_dirs)
    

    lifetimesfile = os.path.join(diagnostic_dir0,agb_mix,set_name,
                                 '_'.join(('tau_cm',name_conv)))
    cm = open(lifetimesfile,'w')
    cm.write('# z mass tauc taum\n')
    for metal_dir in np.array(metal_dirs)[metals]:
        if diagnostic_dir0 != None:
            diagnostic_dir = os.path.join(diagnostic_dir0,agb_mix,set_name,metal_dir)+'/'
            ensure_dir(diagnostic_dir)
            # update kwarg
            kwargs['diagnostic_dir'] = diagnostic_dir
        agb_tracks = get_afile(os.path.join(working_dir,metal_dir),track_identifier)
        agb_tracks.sort()
        ax = plt.axes()
        kwargs['ax'] = ax
        # Paola's formatted tracks
        metallicity = float(metal_dir.split('_')[1].replace('Z',''))
        Y = float(metal_dir.split('_')[-1].replace('Y',''))
        iso_name_conv = '_'.join(('Z%.4f'%metallicity,name_conv))
        isofile = os.path.join(home,isotrack_dir,iso_name_conv)
        isofile_rel_name = os.path.join('isotrack',isotrack_dir,iso_name_conv)
        out = open(isofile,'w')
        out.write('# age(yr) logL logTe m_act mcore c/o period ip Mdot(Msun/yr) logTe X Y CNO dlogTe/dlogL \n')
        print 'found %i tracks'%len(agb_tracks)
        for agb_track in agb_tracks:
            flag = 0
            addpt = []
            # load track
            track = get_numeric_data(agb_track)
            tdata = track.data_array
            ntp = tdata['NTP']
            age = tdata['ageyr']
            phi = tdata['PHI_TP']
            logt = tdata['T_star']
            logl = tdata['L_star']
            dt = tdata['dt']
            co = tdata['CO']
            mdot = tdata['dMdt']
            mass = track.mass
            met = metallicity
            mstar, = np.nonzero((co<=1) & (logl >= 3.3) & (mdot<=-5))
            cstar, = np.nonzero((co>=1) & (mdot<=-5))
            try:
                tauc = np.sum(dt[cstar])/1e6
            except IndexError:
                tauc = 0.
            try:
                taum = np.sum(dt[mstar])/1e6
            except IndexError:
                taum = 0.
            cm.write('%.4f %.3f %.4f %.4f\n'%(met,mass,tauc,taum))
            # get the indices of the thermal pulses.
            if len(ntp) == sum(ntp):
                print 'no tracks!', agb_track
                flag = 1
            if flag == 1: continue
            TPs = get_TP_inds(ntp)
            # The first line in the agb track is 1. This isn't the quiescent...
            phi[0]= -1.
            # The quiescent phase is the the max phase in each TP, i.e., closest to 1.
            
            qs = np.cumsum([np.argmax(phi[TP]) for TP in TPs])
            qs+=range(len(qs))
            #qs = [list(phi).index(np.max(phi[TP])) for TP in TPs]
            Qs = qs
            # when to add an extra point to quiescent track
            if len(Qs) <= 9 and mass < 3.: Qs,addpt = add_points_to_q_track(track,qs)
            Qs = map(int,Qs)
            # find dl/dt of track
            rising,slopes,fits = find_dldt(track,TPs,addpt)
            
            # need to add first line.
            logt0 = logt[0] 
            logl0 = logl[0]
            logt1 = logt[2] 
            logl1 = logl[2]
            slope = (logl1-logl0)/(logt1-logt0)
            Qs.insert(0,0)            
            slopes.insert(0,slope)
            
            # make iso file for trilegal
            make_iso_file(track,Qs,slopes,out)
            # make diagnostic plots and imfr
            if diagnostic_dir != None:
                diag_plots(track,logl,logt,age,slopes,Qs,addpt,rising,fits,**kwargs)
            if kwargs['make_imfr'] == True:
                M_s = AGBTracks.get_col(track,'M_star')
                mfile.write('%f %f %f\n'%(M_s[0],M_s[-1],float(met)))            
        out.close()
        print 'wrote',isofile
        make_local_copy(isofile,dest=make_copy)
        # keep information for tracce file
        isofiles.append(isofile_rel_name)
        Ys.append(Y)
        Zs.append(metallicity)
            
    # make file to link cmd_input to formatted agb tracks
    metfile = make_met_file(tracce_file,Zs,Ys,isofiles)
    print 'wrote',tracce_file
    make_local_copy(tracce_file,dest=make_copy)
    # make cmd_input file that
    cmd_input = write_cmd_input_file(**{'cmd_input_file':cmd_input_file,
                                        'file_tpagb':tracce_file_rel,
                                        'mass_loss':mass_loss})
    print 'wrote',cmd_input_file
    cm.close()
    print 'wrote',lifetimesfile
    plot_cluster_test(lifetimesfile,**kwargs)
    if make_imfr == True:
        mfile.close()
        print 'wrote',ifmr
        plot_ifmr(ifmr)
        
    os.chdir(home)
    return cmd_input_file

if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
    except:
        print do_everything.__doc__
    here = os.getcwd()
    kwargs = load_input(input_file) 
    # Paola's tracks -> trilegal + tests based only on tracks
    cmd_input_file = do_everything(**kwargs)
    
    agb_mix = kwargs['agb_mix']
    set_name = kwargs['set_name']
    track_set = '_'.join((agb_mix,set_name))
    
    diagnostic_dir = kwargs.get('diagnostic_dir')
    # Marco's scripts to run trilegal at age and z
    if kwargs.get('trilegal_diagnostics'): 
        tri_dir = os.path.join(diagnostic_dir,'trilegal_files')
        sfh_dir = os.path.join(tri_dir,'sfh')
        plt_dir = os.path.join(diagnostic_dir,agb_mix,set_name,track_set)    
        trilegal_diagnostics.main(track_set,sfh_dir,tri_dir,plt_dir)
        os.chdir(plt_dir)
        if kwargs.get('google_table'):
            image_location = kwargs.get('image_location')
            GoogleSitesTable.main(image_location)
        os.chdir(here)
    # Scripts to make LF compared to data
    if kwargs.get('IDs'):
        galaxy_tests.run_all(kwargs['IDs'],[track_set+'.dat'])
        pass