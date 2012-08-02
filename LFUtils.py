#
#  LumnFncfromCMD.py
#  
#
#  Created by Philip Rosenfield on 11/22/10.
#
import sys,os,re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib import rc
rc('text',usetex=True)

from scipy.stats import ks_2samp
# Phil's functions:
from MatchUtils import readbintab
from GenUtils import Table,read_table,get_afile,discrete_colors
from Astronomy import get_tab5_trgb_Av_dmod
from TPAGBparams import *
from TrilegalUtils import get_stage_label, get_label_stage
from matplotlib.patheffects import withStroke

import brewer2mpl
myeffect = withStroke(foreground="w", linewidth=3)
kwargs = dict(path_effects=[myeffect])

def read_irtable_data(tab=os.path.join(table_src,'IR_NAGBs.dat')):
    '''OLD'''
    lines = open(tab,'r').readlines()
    col_keys = lines[0].replace('#','').strip().split()
    d = {}
    for line in lines:
        if line.startswith('#'): continue
        ID = line.split('_')[1]
        if ID not in d.keys(): d[ID] = {}
        data = line.strip().split()[1:]
        for i,v in enumerate(data):
            j = i+1
            if col_keys[j] not in d[ID].keys(): d[ID][col_keys[j]] = {}
            d[ID][col_keys[j]] = GenUtils.is_numeric(v)
    return d

def get_irtable_data(Target):
    tab = read_irtable_data()
    return tab[Target]

def read_mettable():
    tab = '/Users/phil/research/TP-AGBcalib/SNAP/tables/IR_NAGBs.dat'
    dtype=[('fitstable', '|S46'), ('IR_TRGB', '<f8'), ('N_AGB', '<f8'), ('logOH', '<f8'), ('OHerr', '<f8'), ('Z', '<f8'), ('FeH', '<f8'), ('Ttype', '<f8')]
    table = np.genfromtxt(tab,dtype=dtype,delimiter=',')
    return table

def get_key_fromtable(ID,key):
    tab = read_mettable()
    names = list(tab['fitstable'])
    i, = [names.index(i) for i in names if ID in i]
    return tab[key][i]


def get_trgb_ir_nAGB(Target):
    trgb_ir = get_key_fromtable(Target,'IR_TRGB')
    nAGB = get_key_fromtable(Target,'N_AGB')
    '''
    table = open(os.path.join(table_src,'IR_NAGBs.dat'),'r')
    i = 0
    for line in table.readlines():
        if line.startswith('#'): continue
        a = line.strip().split()
        fitsfile,strgb_ir,snAGB = a[0:3]
        ID = fitsfile.split('_')[1]
        if Target != ID: continue
        i = 1
        trgb_ir = float(strgb_ir)
        nAGB = int(snAGB)
    if i == 0: 
        print 'No listing of',Target,'in',table_src
        trgb_ir = 99.99
        nAGB = 0
    '''
    return trgb_ir,nAGB

def brighter_trgb(mag2,trgb,inds=None):
    # number of stars brighter than trgb, make sure mag2 is 
    # the same filter as trgb!
    if inds==None:
        i, = np.nonzero(mag2 < trgb)
    else:
        j, = np.nonzero(mag2 < trgb)
        i = list(set(j) & set(inds))
    return i,len(i)

def get_chi2(obs,exp):
    return np.sum((obs-exp)**2/(exp))

def nbetween(arr,mdim,mbrt):
    n, = np.nonzero((arr<mdim) & (arr>mbrt))
    return n.size

def calc_LF(gal,sgal,maglims,res=0.1):
    '''
    usage: galaxy object, simgalaxy object, offset (mags below trgb), LF bin resolution
    
    returns:
    ind: indices of simulation cmd that are picked so they mimic the number of rgb stars
    between trgb and trgb+offset

    p_value: ks test probability comparing simulated stars brighter than trgb to data of
    the same

    normalization: scaling used to determine ind and scale LFs (added attributes to 
    objects)
    '''
            
    mag = gal.mag2
    dim,bright = maglims
    # when hand picking cmd_regions, RGB includes AGB!
    gal.irgb = gal.stage_inds('RGB')
        
    smag = sgal.ast_mag2
    sgal.irgb = sgal.stage_inds('RGB')
    
    gal.iAGB,gal.nAGB = brighter_trgb(mag,bright,inds=gal.irgb)
    
    #res = 0.1 # LF resolution in dex
    nbins = (mag.max() - mag.min()) / res
    
    # The RGB stars used for normalization
    gal.nRGB = float(nbetween(mag[gal.irgb],dim,bright))
    #nRGB = nbetween(mag,norm,trgb)
    
    # The Simulated RGB stars
    sgal.nRGB = float(nbetween(smag[sgal.irgb],dim,bright))
    
    normalization = gal.nRGB/sgal.nRGB # here is the normalization!!
    if normalization > 0.75:
        print normalization
        print 'Too few model stars'
    
    # the LFs
    gal.LF,gal.bins = np.histogram(mag,nbins)
    sgal.LF,sgal.bins = np.histogram(smag,bins=gal.bins)
    
    # normalize the simulated LF
    sgal.LFn = sgal.LF*normalization 
    
    # for plotting random simulated stars to match number of obs.
    # jesus fucking christ. should this be smag.size or smag[sgal.irgb].size?!
    #rands = np.random.random(smag.size)
    rands = np.random.random(smag[sgal.irgb].size)
    ind, = np.nonzero(rands<normalization)
    
    # number of simulated RGB stars (now normalized to match the relative number of RGBs)
    sgal.nRGBn = nbetween(smag[ind],dim,bright)

    # the indices and number of simulated (and normalized) stars brighter than the trgb
    sgal.iAGBn,sgal.nAGBn = brighter_trgb(smag[ind],bright)
    
    KS_D,p_value = ks_2samp(mag[gal.iAGB],smag[sgal.iAGBn])
    #print ps_nNorm,nNorm,s_nNorm,np.sqrt(nNorm),ps_nNorm-nNorm
    #print ps_nB_AGBm,nB_AGB
    # want to add some lines for the flux and the mass loss rates.
    return ind,p_value,normalization

    
def setup_lfplot(gal):
            
    # plot limits
    left, width = 0.1, 0.312
    bottom, height = 0.1, 0.8
    left_m = left+width+0.01 # for model
    left_h = left+2*width+0.02 # for hist (width is set hist_axis)
    
    # plot and fig sizes
    fig = plt.figure(1, figsize=(8,8))
    
    data_axis = [left, bottom, width, height]
    model_axis =[left_m, bottom, width, height] 
    hist_axis = [left_h, bottom, 0.2, height]
    
    axData = plt.axes(data_axis)
    axModel = plt.axes(model_axis)    
    axHist = plt.axes(hist_axis)
    
    # titles
    axData.set_title(r'$\rm{data}$',color = 'black')
    axModel.set_title(r'$\rm{model}$',color = 'red')
    axData.set_xlabel(r'$\rm{%s-%s}$'%(gal.filter1,gal.filter2),size=20)
    axData.set_ylabel(r'$\rm{%s}$'%gal.filter2, size=20)
    axModel.set_xlabel(axData.get_xlabel(), size=20)
    axHist.set_xlabel( r'$\#$',size=20)
    
    for ax in [axData,axModel,axHist]:
        ax.set_xlim( (-0.5, gal.color.max()) ) # model and data x limits here
        ax.set_ylim( (25.5, 18.5) ) # set all y limits here
    axHist.set_xlim( (1, 50000) )
    
    return fig,axData,axModel,axHist

def plot_lines(axs,xrange,yval):
    x=[ax.plot((xrange),(yval,yval),color='black') for ax in axs]
    return
    
def plot_numbs(ax,item,xpos,ypos,**kwargs):
    x= ax.annotate(r'$%i$' % item,xy=(xpos,ypos),ha='left',size=20,**kwargs)
    return

def diagnostic_cmd(sgal,trgb,figname=None):
    ustage = np.unique(sgal.stage)
    nplots = ustage.size+1.
    cols = brewer2mpl.get_map('RdBu','diverging',len(ustage)).mpl_colors
    subplots_kwargs = {'sharex':1,'sharey':1,'figsize':(12, 8)}
    fig, axs = setup_multiplot(nplots,**subplots_kwargs)
    
    for ax in axs.ravel():
        ax.set_xlim(-0.5,sgal.color.max())
        ax.set_ylim(25.5,18.5)
                
    j=0
    for color,mag2 in zip((sgal.color,sgal.ast_color),(sgal.mag2,sgal.ast_mag2)):
        ax0,cols = colorplot_by_stage(axs.ravel()[0],color,mag2,'.',sgal.stage,cols=cols)
        i=0
        for ax,st in zip(axs.ravel()[1:],ustage): 
            plot_lines([ax],ax.get_xlim(),trgb)
            label = get_label_stage(int(st))
            ind = sgal.stage_inds(label)
            if len(ind)==0: continue
            ax.plot(color[ind],mag2[ind],'.',color=cols[i],mew=0,label=label)
            ax.set_title(label,**{'color':cols[i]})
            i+=1
        if figname:
            if j==0: 
                extra = ''
            else:
                extra = '_spread'
            plt.savefig(figname.replace('.png','%s.png'%extra))
            print 'wrote %s'%figname.replace('.png','%s.png'%extra)

        else:
            plt.show()
        j+=1
    plt.close()
    return
    
def setup_multiplot(nplots,**subplots_kwargs):
    nx = np.round(np.sqrt(nplots))
    nextra = nplots-nx**2
    ny = nx
    if nextra > 0: ny += 1
    nx = int(nx)
    ny = int(ny)
    
    fig, axs = plt.subplots(nx,ny,**subplots_kwargs)
    
    return fig, axs
    
def colorplot_by_stage(ax,x,y,marker,stages,cols=None):
    # inds from calc_LFIR are based on only resolved stars.
    
    if cols == None:
        cols = discrete_colors(len(np.unique(stages)))
    for i,s in enumerate(np.unique(stages)):
        ind, = np.nonzero(stages == s)
        if ind.size == 0: continue
        ax.plot(x[ind],y[ind],marker,color=cols[i],mew=0,label=get_label_stage(int(s)))
    plt.legend()
    return ax,cols

def make_title(gal,fig):
    text_kwargs = {'ha':'center','va':'top','size':20}
    if np.isfinite(gal.z):
        title =r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=%.4f$'%(gal.target,gal.trgb,gal.z)
    else:
        title = r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=...$'%(gal.target,gal.trgb)
        
    x = fig.text(0.5,0.96,title,**text_kwargs)
    return
    
def plot_LFIR(gal,sgal,ind,p_value,maglims):
    dim,bright = maglims

    mstars = list(set(sgal.imstar) & set(ind))
    cstars = list(set(sgal.icstar) & set(ind))
    nRGBs = len(sgal.ast_color[ind]) - len(mstars) - len(cstars)
   
    mhist,b = np.histogram(sgal.mag2[mstars],bins=sgal.bins)
    chist,b = np.histogram(sgal.mag2[cstars],bins=sgal.bins)
    
    fig,axData,axModel,axHist = setup_lfplot(gal)    

    make_title(gal,fig)
    
    # plot data
    axData.plot(gal.color,gal.mag2,'.',mew=0,color='grey',mec='grey')
    axData.plot(gal.color[gal.irgb],gal.mag2[gal.irgb],'.',mew=0,color='black')
    
    axModel.plot(sgal.ast_color[ind],sgal.ast_mag2[ind],'.',
                 color='red',mew=0,label=r'$RGB=%i$'%(nRGBs))
    axModel.plot(sgal.ast_color[mstars],sgal.ast_mag2[mstars],'.',
                 color='darkgreen',mew=0,label=r'$M=%i$'%(len(mstars)))
    axModel.plot(sgal.ast_color[cstars],sgal.ast_mag2[cstars],'.',
                 color='darkblue',mew=0,label=r'$C=%i$'%(len(cstars)))

      
      
    axModel.legend(frameon=False,loc=2,numpoints=1)
    
    [plot_lines([axData,axModel],axData.get_xlim(),m) for m in (dim,bright)]
    # also?
    # plot_lines([axData,axModel],gal.color.min(),gal.color.max(),trgb)
    
    kwargs['color'] = 'black'
    text_offset =0.02
    xpos = axData.get_xlim()[0]+2*text_offset
    yposs = np.asarray(maglims)-text_offset
    plot_numbs(axData,gal.nAGB,xpos,yposs[1],**kwargs)
    plot_numbs(axData,gal.nRGB,xpos,yposs[0],**kwargs)
    
    kwargs['color'] = 'red'
    xpos = axModel.get_xlim()[0]+text_offset
    plot_numbs(axModel,sgal.nAGBn,xpos,yposs[1],**kwargs)
    plot_numbs(axModel,sgal.nRGBn,xpos,yposs[0],**kwargs)
    
    # plot hists
    #hist,bins,patches = axHist.hist(mag2,nbins,histtype='step',orientation='horizontal',log=True,color='black')
    #s_hist,s_bins,s_patches = axHist.hist(s_mag2,s_nbins,histtype='step',orientation='horizontal',log=True,color='red')
    #axHist.semilogx(hist,bins[:-1],drawstyle='steps',color='black')

    
    axHist.semilogx(mhist,b[:-1],drawstyle='steps',color='darkgreen')
    axHist.semilogx(chist,b[:-1],drawstyle='steps',color='darkblue')
    axHist.semilogx(gal.LF,b[:-1],drawstyle='steps',color='black')
    axHist.semilogx(sgal.LFn,b[:-1],drawstyle='steps',color='red')
    kwargs['color']='black'
    axHist.annotate(r'$p=%.3f$' % p_value,xy=(.9,.95), xycoords='axes fraction',ha='right',size=20,**kwargs)
    
    # no formatters on mid and right plots
    nullfmt   =  NullFormatter() # no labels
    [ax.yaxis.set_major_formatter(nullfmt) for ax in (axHist,axModel)]
    
    filename= '_'.join((gal.target,sgal.model,gal.filter1,gal.filter2))
    #plt.savefig(os.path.join(plt_dir,filename+'_LF.ps'))
    plt.savefig(os.path.join(plt_dir,filename+'_LF.png'))
    print 'Wrote ',filename+'_LF.png'
    plt.close()
    return
    
def write_cmd_forJason(filename,color,mag,filt1,filt2):
    # For Jason
    if not os.path.isdir('JASON'): os.mkdir('JASON')
    outfile = open('JASON/'+filename+'_cmd.dat','w')
    outfile.write('#%s-%s %s\n'%(filt1,filt2,filt2))
    [outfile.write('%f %f\n'%(color[i],mag[i])) for i in range(len(mag))]
    outfile.close()
    print 'wrote JASON/'+filename+'_cmd.dat'
    return


def make_galaxy_params_table(Target):
    for line in open('/Users/Phil/research/Italy/WFC3SNAP/PHIL/table.dat','r').readlines():
        if line.startswith('#'): continue
        Target,dmod,av,camera,mag50,filt50 =  line.split()
        trgb,av,dmod = get_tab5_trgb_Av_dmod(Target)
        if Target == 'NGC0404-DEEP': Target = 'NGC404'
        trgb_ir,Nagb = get_trgb_ir_nAGB(Target)
        print Target,dmod,av,camera,mag50,filt50,trgb,trgb_ir
        
if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print 'usage: (opt or ir)'
        print 'python LFUtils.py band Target trilegal_output_AST_corrected model'
        print 'python LFUtils.py ir IC2574-SGS output_IC2574-SGS_spread_IR.dat cmd_input_gi10_rev.dat'
    band = sys.argv[1]
    Target = sys.argv[2]
    trilegal_out = sys.argv[3]
    model  = sys.argv[4]
    
    if band == 'opt':
        fits = get_afile(fits_src,'*'+'*'.join((Target,'trim','.fits')))[0]
        trgb = get_tab5_trgb_Av_dmod(Target)[0]
    elif band == 'ir':
        # read data
        fits = get_afile(fits_src,'*'+'*'.join((Target,'IR','.fits')))[0]
        trgb,nAGB = get_trgb_ir_nAGB(Target)
    else: 
        print 'choose opt or ir'
        sys.exit()
    if trilegal_out == None:
        trilegal_out = get_afile(output_src,'*'+'*'.join(('ast',Target,model))+'*')[0] # and ast...
        print 'no trilegal output file supplied, using',trilegal_out
    
    ## INPUT ##
    # read data
    f = readbintab(fits)
    filt1 = f['names'][-1].split('.')[0].split('_')[-2]
    filt2 = f['names'][-1].split('.')[0].split('_')[-1]
    print 'fitstable =',f['names'][-1]
    mag1 = f['Mag1']
    mag2 = f['Mag2']
    color = mag1-mag2
    
    # read model
    print 'reading', trilegal_out
    synthcmd = get_numeric_data(trilegal_out)
    print 'ok.'
    s_mag1 = synthcmd.get_col(filt1) + synthcmd.get_col('diff_'+filt1.strip())
    s_mag2 = synthcmd.get_col(filt2) + synthcmd.get_col('diff_'+filt2.strip())
    
    s_mag1 = s_mag1[np.nonzero(abs(synthcmd.get_col('diff_'+filt1.strip())) < 90.)[0]]
    s_mag2 = s_mag2[np.nonzero(abs(synthcmd.get_col('diff_'+filt1.strip())) < 90.)[0]]
    s_color = s_mag1-s_mag2
    
    Norm = trgb + 1.5
    ind,nB_AGB,nNorm,ps_nNorm,ps_nB_AGBm,hist,bins,s_hist_normed,p_value,normalization = calc_LF(mag2,s_mag2,Norm,trgb)
    filename= '_'.join((Target,model,filt1,filt2))
    write_cmd_forJason(filename,s_color[ind],s_mag2[ind],filt1,filt2)
    plot_LFIR(Target,model,trgb,Norm,filt1,filt2,color,mag2,s_color,s_mag2,ind,nB_AGB,nNorm,ps_nNorm,ps_nB_AGBm,hist,bins,s_hist_normed,p_value,normalization)
    
    

