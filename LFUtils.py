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
import GenUtils
from Astronomy import get_tab5_trgb_Av_dmod
from TPAGBparams import *
from TrilegalUtils import get_stage_label, get_label_stage
from matplotlib.patheffects import withStroke

import brewer2mpl
myeffect = withStroke(foreground="w", linewidth=3)
kwargs = dict(path_effects=[myeffect])

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
    return trgb_ir,nAGB

def brighter(mag2,trgb,inds=None):
    # number of stars brighter than trgb, make sure mag2 is 
    # the same filter as trgb!
    i, = np.nonzero(mag2 < trgb)
    if inds!=None: i = np.intersect1d(i,inds)
    return i

def get_chi2(obs,exp):
    return np.sum((obs-exp)**2/(exp))

def between(arr,mdim,mbrt,inds=None):
    i, = np.nonzero((arr<mdim) & (arr>mbrt))
    if inds!=None:
        i = np.intersect1d(i,inds)
    return i

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
    
    # data
    mag = gal.mag2
    color = gal.color
    dim,bright = maglims
    # when hand picking cmd_regions, RGB includes AGB!
    gal.irgb = gal.stage_inds('RGB')
    gal.iagb = brighter(mag,gal.trgb,inds=gal.irgb)
    
    # simulation
    smag = sgal.ast_mag2
    scolor = sgal.ast_color
    sgal.irgb = sgal.stage_inds('RGB')
        
    # The RGB stars used for normalization
    gal.rgb_norm = between(mag,dim,bright,inds=gal.irgb)
    
    # The Simulated RGB stars
    #!! Not used!!
    #sgal.rgb_in_norm = between(smag,dim,bright,inds=sgal.irgb)
    
    # all sim stars in cmd region of rgb data norm stars.
    sgal.norm = GenUtils.inside(color[gal.rgb_norm],mag[gal.rgb_norm],
                                scolor,smag,)
    
    # here is the normalization!!
    sgal.normalization = float(gal.rgb_norm.size)/float(sgal.norm.size)
    if sgal.normalization > 0.75:
        print sgal.normalization
        print 'Too few model stars'
    
    # LF resolution in dex
    nbins = (mag.max() - mag.min()) / res
    
    # the LFs
    gal.LF,gal.bins = np.histogram(mag,nbins)
    sgal.LF,sgal.bins = np.histogram(smag,bins=gal.bins)
    
    # normalize the simulated LF
    sgal.LFn = sgal.LF*sgal.normalization 
    
    # for plotting random simulated stars to match number of obs.
    rands = np.random.random(smag.size)
    ind, = np.nonzero(rands < sgal.normalization)
    sgal.rel_ind = ind
    
    # simulated normalized RGB stars 
    sgal.rel_rgb = between(smag[ind],dim,gal.trgb)
    
    # simulated normalized stars brighter than trgb.
    sgal.rel_agb = brighter(smag[ind],gal.trgb)
    
    KS_D,p_value = ks_2samp(mag[gal.iagb],smag[sgal.rel_agb])

    # want to add some lines for the flux and the mass loss rates.
    return p_value

    
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
    fig, (axs) = setup_multiplot(nplots,**subplots_kwargs)
    
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
            ax.plot(color[ind],mag2[ind],'.',color=cols[i],mew=0)
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
    
def plot_LFIR(gal,sgal,p_value,maglims):
    dim,bright = maglims

    mstars = list(set(sgal.imstar) & set(sgal.rel_ind))
    cstars = list(set(sgal.icstar) & set(sgal.rel_ind))
    #nRGBs = len(sgal.ast_color[sgal.rel_agb]) - len(mstars) - len(cstars)
   
    mhist,b = np.histogram(sgal.mag2[mstars],bins=sgal.bins)
    chist,b = np.histogram(sgal.mag2[cstars],bins=sgal.bins)
    
    fig,axData,axModel,axHist = setup_lfplot(gal)    

    make_title(gal,fig)
    
    # plot data
    axData.plot(gal.color,gal.mag2,'.',mew=0,color='grey',mec='grey')
    # data used for normalization
    axData.plot(gal.color[gal.rgb_norm],gal.mag2[gal.rgb_norm],'.',mew=0,color='black')
    
    # plot model
    axModel.plot(sgal.ast_color[sgal.rel_ind],sgal.ast_mag2[sgal.rel_ind],'.',
                 color='red',mew=0)#,label=r'$RGB=%i$'%(nRGBs))                 
                 
    # model used for normalization
    norm_inds = list(set(sgal.norm) & set(sgal.rel_ind))
    axModel.plot(sgal.ast_color[norm_inds],sgal.ast_mag2[norm_inds],'.',
                 color='black',mew=0)#,label=r'$RGB=%i$'%(len(norm_inds)))

    axModel.plot(sgal.ast_color[mstars],sgal.ast_mag2[mstars],'.',
                 color='darkgreen',mew=0,label=r'$M=%i$'%(len(mstars)))

    axModel.plot(sgal.ast_color[cstars],sgal.ast_mag2[cstars],'.',
                 color='darkblue',mew=0,label=r'$C=%i$'%(len(cstars)))
    
    axModel.legend(frameon=False,loc=2,numpoints=1)
    
    #[plot_lines([axData,axModel],axData.get_xlim(),m) for m in (dim,bright)]
    # also?
    plot_lines([axData,axModel],axData.get_xlim(),gal.trgb)
    
    kwargs['color'] = 'black'
    text_offset = 0.02
    xpos = axData.get_xlim()[0]+2*text_offset
    yposs = np.asarray(maglims)-text_offset
    ypos = gal.trgb-text_offset
    plot_numbs(axData,gal.iagb.size,xpos,ypos,**kwargs)
    plot_numbs(axData,gal.rgb_norm.size,xpos,yposs[0],**kwargs)
    
    kwargs['color'] = 'red'
    xpos = axModel.get_xlim()[0]+text_offset
    plot_numbs(axModel,sgal.rel_agb.size,xpos,ypos,**kwargs)
    plot_numbs(axModel,len(norm_inds),xpos,yposs[0],**kwargs)
    
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
    
    

