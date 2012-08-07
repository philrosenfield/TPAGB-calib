import operator
import cmdUtils
from TrilegalUtils import get_stage_label
import LFUtils
import mk_sims
from multiprocessing import Pool
import itertools
import time
from TPAGBparams import *
import GenUtils
import numpy as np
import pdb; pdb.set_trace()
import matplotlib.pyplot as plt
import GoogleSitesTable as gst

def tag_cmds(IDs):
    '''
    only need to do this once, hopefully. Just went though and added integer tags to each
    fits table so that it would be like trilegal output with the -l option.
    '''
    seqs = ['MS','RHeB','RGB']
    if type(IDs) == str:
        fits = GenUtils.get_afile(fits_src,'*'+'*'.join((IDs,'IR','.fits')))
    else:
        fits = np.squeeze([GenUtils.get_afile(fits_src,'*'+'*'.join((ID,'IR','.fits'))) for ID in IDs]) 
    cmdUtils.define_color_mag_region(fits,seqs)
    return
    
def run_all(IDs,models):
    pool=Pool()
    
    res=[]
    for ID,model in itertools.product(IDs,models):
        res.append(pool.apply_async(main,(ID,model),))

    for r in res:
        r.get()
    
    return
    
def read_tagged_phot(tagged_file):
    '''
    reads a file created by cmdUtils.define_color_mag_region. ascii with 7 columns.
    '''
    if type(tagged_file) == str:
        #print 'reading %s'%tagged_file
        cols = ['ra','dec','mag1','mag2','mag1err','mag2err','stage']
        fits = np.genfromtxt(tagged_file,names=cols)
    else:
        fits = tagged_file
    return fits
    
def get_trgb_fitsname(ID,band):
    if band == 'opt':
        fits = GenUtils.get_afile(fits_src,'*'+'*'.join((ID,'trim','.fits')))[0]
        trgb = get_tab5_trgb_Av_dmod(ID)[0]
    elif band == 'ir':
        # read data
        fits, = GenUtils.get_afile(fits_src,'*'+'*'.join((ID,'IR','.fits')))
        trgb = LFUtils.get_trgb_ir_nAGB(ID)[0]
    else: 
        print 'choose opt or ir'
        sys.exit()
    return trgb,fits

'''
def initial_object_masses(ID):
    d = {"SCL-DE1": 5e+06,
            "NGC2403-HALO-6": ,
            "NGC7793-HALO-6": ,
            "UGC-04459": ,
            "UGC-4305-1": ,
            "UGC-4305-2": ,
            "UGC-5139": ,
            "DDO78": ,
            "DDO82": ,
            "KDG73": ,
            "KKH37": ,
            "NGC0300-WIDE1": ,
            "NGC404": ,
            "NGC2403-DEEP": ,
            "NGC2976-DEEP": ,
            "NGC3741": ,
            "NGC4163": ,
            "UGC8508": ,
            "UGCA292": ,
            "NGC3077-PHOENIX": ,
            "IC2574-SGS": ,
            "HS117": ,
            "DDO71":
            }
    return d[ID]
'''

class galaxies(object):
    '''
    THIS IS FROM BR RATIO CODE. I'M NOT SURE IT'S GENERAL ENOUGH TO PUT INTO
    ~/research/python SO I COPIED IT...
    
    I made summary... 
    wrapper for lists of galaxy objects, each method returns lists, unless they 
    are setting attributes.
    '''
    def __init__(self, galaxy_objects):
        self.galaxies = galaxy_objects
        galaxies.summary(self,'target','filter1','filter2')
        # maybe I can avoid simgalaxies with self.__class__?
        #self.zs = np.unique([np.round(g.z, 4) for g in self.galaxies])
        #self.photsyss =  np.unique(g.photsys for g in self.galaxies)
    
    def summary(self,*attrs):
        '''
        similar to squish, but has np.unique.
        '''
        for attr in attrs:
            new_list = [g.__getattribute__(attr) for g in self.galaxies]
            new_val = np.unique(new_list)
            self.__setattr__('%ss'%attr,new_val)
        
    def all_stages(self, *stages):
        '''
        adds the indices of any stage as attributes to galaxy.
        If the stage isn't found, -1 is returned.
        '''
        [g.all_stages(*stages) for g in self.galaxies]
        return
    
    def squish(self, *attrs):
        '''
        concatenates an attribute or many attributes and adds them to galaxies 
        instance.
        No slicing, so not sure how it will be useful besides Color Mag2.
        '''
        for attr in attrs:
            new_list = [g.__getattribute__(attr) for g in self.galaxies]
            new_val = np.concatenate(new_list)
            self.__setattr__('%ss'%attr, new_val)
            
    def finite_key(self, key):
        return [g for g in self.galaxies if np.isfinite(g.__dict__[key])]
    
    def select_on_key(self, key, val):
        ''' ex filter2 == F814W works great with strings or exact g.key==val.
        rounds z to four places, no error handling.
        '''
        key = key.lower()
        if key == 'z':
            gs = [g for g in self.galaxies if np.round(g.__dict__[key], 4)==val]
        else:
            gs = [g for g in self.galaxies if g.__dict__[key]==val]
        return gs
            
    def group_by_z(self):
        zsf = self.zs[np.isfinite(self.zs)]
        zsn = self.zs[np.isnan(self.zs)]
        d = {}
        for z in zsf:
            key = 'Z%.4f'%z
            d[key] = galaxies.select_on_key(self, 'z', z)
        
        d['no z'] = [g for g in gals if np.isnan(g.z)]
        return d

    def intersection(self, **kwargs):
        '''
        ex kwargs = {'filter2':'F814W', 'filter1':'F555W'}
        will return a list of galaxy objects that match all kwarg values.
        '''
        gs_tmp = self.galaxies
        gs = [galaxies.select_on_key(self, k, v) for k, v in kwargs.items()]
        for i in range( len(gs) ):
            gs_tmp = list(set(gs_tmp) & set(gs[i]))
        return gs_tmp
    
    def __str__(self):
        for g in self.galaxies:
            print g.__str__()
        return ''


class galaxy(object):
    def __init__(self, ID, band):
        self.name = get_trgb_fitsname(ID,band)[1]
        self.data = read_tagged_phot(cmdUtils.GenUtils.replace_ext(self.name,'.dat'))
        self.base = os.path.split(self.name)[0]
        self.target = cmdUtils.get_fileinfo(self.name)[0]
        self.filter1 = cmdUtils.get_fileinfo(self.name)[1]
        self.filter2 = cmdUtils.get_fileinfo(self.name)[2]
        self.mag1 = self.data['mag1']
        self.mag2 = self.data['mag2']
        self.color = self.data['mag1']-self.data['mag2']
        self.stage = self.data['stage']
        self.trgb = get_trgb_fitsname(ID,band)[0]
        self.z = LFUtils.get_key_fromtable(ID,'Z')
        
    def stage_inds(self,stage_name):
        return get_stage_inds(self.data,stage_name)
    
    def delete_data(self):
        data_names = ['data','mag1','mag2','color','stage']
        [self.__delattr__(data_name) for data_name in data_names]
        
    def all_stages(self, *stages):
        '''
        adds the indices of some stage as an attribute.
        '''
        for stage in stages:
            i = stage_inds(self.stage, stage)
            self.__setattr__('i%s'%stage.lower(), i)
        return
                
class simgalaxy(object):
    def __init__(self,trilegal_out,filter1,filter2):
        self.base,self.name = os.path.split(trilegal_out)
        self.data = GenUtils.read_table(trilegal_out)
        self.filter1 = filter1
        self.filter2 = filter2
        self.target = self.name.split('_')[2]
        self.mag1 = self.data.get_col(self.filter1)
        self.mag2 = self.data.get_col(self.filter2)
        self.stage = self.data.get_col('stage')
        
        simgalaxy.load_ast_corrections(self)
        
        data_to_slice = ['mag1','mag2','stage','ast_mag1','ast_mag2']
        slice_inds = self.rec
        simgalaxy.slice_data(self,data_to_slice,slice_inds)
            
        self.ast_color = self.ast_mag1-self.ast_mag2
        self.color = self.mag1-self.mag2
        simgalaxy.load_ic_mstar(self)
        
    def load_ast_corrections(self):
        diff1 = self.data.get_col('diff_'+self.filter1)
        diff2 = self.data.get_col('diff_'+self.filter2)
        recovered1, = np.nonzero(abs(diff1)<90.)
        recovered2, = np.nonzero(abs(diff2)<90.)
        self.rec = list(set(recovered1) & set(recovered2))
        self.ast_mag1 = self.mag1+diff1
        self.ast_mag2 = self.mag2+diff2
    
    def slice_data(self,data_to_slice,slice_inds):
        '''
        slice already set attributes by some index list.
        '''
        [self.__setattr__(d,self.__dict__[d][slice_inds]) for d in data_to_slice]
        
    def mix_modelname(self,model):
        self.mix,self.model_name = get_mix_modelname(model)
    
    def delete_data(self):
        '''
        for wrapper functions, I don't want gigs of data stored when they
        are no longer needed.
        '''
        data_names = ['data','mag1','mag2','color','stage','ast_mag1',
                      'ast_mag2','ast_color','rec']
        [self.__delattr__(data_name) for data_name in data_names]

    def stage_inds(self,stage_name):
        return np.nonzero(self.stage == get_stage_label(stage_name))[0]

    def load_ic_mstar(self):
        co = self.data.get_col('C/O')[self.rec]
        lage = self.data.get_col('logAge')[self.rec]
        mdot = self.data.get_col('logML')[self.rec]
        logl = self.data.get_col('logL')[self.rec]
        
        self.imstar, = np.nonzero((co<=1) & 
                                  (logl >= 3.3) & 
                                  (mdot<=-5) & 
                                  (self.stage == get_stage_label('TPAGB')))
                                  
        self.icstar, = np.nonzero((co>=1) & 
                                  (mdot<=-5) & 
                                  (self.stage == get_stage_label('TPAGB')))
    
    def all_stages(self, *stages):
        '''
        adds the indices of some stage as an attribute.
        '''
        for stage in stages:
            i = stage_inds(self.stage, stage)
            self.__setattr__('i%s'%stage.lower(), i)
        return

def get_mix_modelname(model):
    mix = model.split('.')[0].split('_')[2]
    model_name = '_'.join(model.split('.')[0].split('_')[3:])
    return mix,model_name

def load_galaxy_tagged(ID,band):
    trgb,fitsname = get_trgb_fitsname(ID,band) 
    tagged_fits = read_tagged_phot(GenUtils.replace_ext(fitsname,'.dat'))
       
    target, filt1,filt2 = cmdUtils.get_fileinfo(fitsname)
    mag1 = tagged_fits['mag1']
    mag2 = tagged_fits['mag2']
    stage = tagged_fits['stage']
    
    color = mag1-mag2
    
    return tagged_fits
    
def get_stage_inds(fits,stage_name):
    inds, = np.nonzero(fits['stage'] == get_stage_label(stage_name))
    return inds

def make_normalized_simulation(ID,model,**kwargs):
    '''
    Wrapper for running trilegal and LFUtils.calc_LF 
    Will continue to run trilegal until galaxy is high enough mass to have 
    proper normalization.
    
    many attributes of galaxy and simgalaxy instance are set in LFUtils.calc_LF.
    
    inputs:
    ID of target: string
    cmd_input model to run trilegal: string
    
    optional params:
    over_write: bool [False]
        force running of trilegal (will run if no output file exists)
    maglims: list, tuple or string. ['trgb']
        either the dim,bright maglims or 'trgb'
    offsets: list, tuple: (dim,bright) default: [(1.5,0.5)]
        must exist if maglims='trgb'
        
        maglims will be trgb+dim, trgb+bright
    object_mass: float [5e6]
        initial mass to run trilegal. 
        
    returns:
        gal instance
        sgal instance
        ks test p_value
        maglims
    '''
    # To be compatible with Jason, might work with 'ir' or 'opt' should only
    # call different file finding routines.
    band = kwargs.get('band','ir')
    # load kwargs
    over_write = kwargs.get('over_write',False)
    object_mass = kwargs.get('object_mass',5e6)
    maglims = kwargs.get('maglims','trgb')
    offsets = kwargs.get('offsets',(1.5,0.5))
    
    # data galaxy instance
    gal = galaxy(ID,band)
    # set maglim if not already set.
    if maglims == 'trgb':
        maglims = (gal.trgb+offsets[0],gal.trgb+offsets[1])

    # initializations
    go = 0
    normalization = 1e9
    # arbitrary normalization factor set to 0.75... 1.0 is upper limit, 
    # didn't think it's necessary to make this a kwarg...
    while normalization >.75:
        # object_mass increase factor is hard coded, could be kwarg.
        if go != 0: 
            object_mass = object_mass*5.
            # reset over_write for while loop
            over_write = 1
        go +=1
        
        # run trilegal
        if over_write:
            print 'Trying %s %s, Mass = %g, Attempt %i'%(ID,model,object_mass,go)
        trilegal_out = mk_sims.mk_sims(ID,model,
                                       object_mass=object_mass,
                                       over_write=over_write)
        
        # load simgalaxy
        sgal = simgalaxy(trilegal_out,gal.filter1,gal.filter2)
        sgal.ID = ID
        sgal.model = model
        sgal.mix_modelname(model)
        
        # calcuate the LFs
        p_value = LFUtils.calc_LF(gal,sgal,maglims)
        
        # update normalization
        normalization = sgal.normalization
        print 'normalization',normalization
    
    if over_write: 
        sgal.object_mass = object_mass
        print '\n %s necessary object_mass = %g\n'%(gal.name,object_mass)
        
    return gal,sgal,p_value,maglims 


def all_IDs():
    IDs=["SCL-DE1",
        "NGC2403-HALO-6",
        "NGC7793-HALO-6",
        "UGC-04459",
        "UGC-4305-1",
        "UGC-4305-2",
        "UGC-5139",
        "DDO78",
        "DDO82",
        "KDG73",
        "KKH37",
        "NGC0300-WIDE1",
        "NGC404",
        "NGC2403-DEEP",
        "NGC2976-DEEP",
        "NGC3741",
        "NGC4163",
        "UGC8508",
        "UGCA292",
        "NGC3077-PHOENIX",
        "IC2574-SGS",
        "HS117",
        "DDO71"]
    return IDs

def main(make_plots=False,publish_plots=False):

    # load all targets
    IDs = all_IDs()
    #IDs = [IDs[0]]
    
    # load all models
    models = ["cmd_input_CAF09_S_SCS.dat",
              "cmd_input_CAF09_S_SCSFG.dat",
              "cmd_input_CAF09_S_SCSFG_ETA2.dat"]
    
    outfile = 'result_tab.dat'
    if os.path.isfile(outfile): 
        print outfile,'exists. appending.'
        out = open(outfile,'a')
    else:
        out = open(outfile,'w')
        out.write('# ID model p_value NRGB_data NAGB_data NRGB_model NAGB_model  \n')
    
    for ID,model in itertools.product(IDs,models):
        gal, sgal, p_value, maglims = make_normalized_simulation(ID,model)
        write_spread_catalog(sgal)
        
        if make_plots:
            fig_loc = os.path.join(plt_dir,sgal.mix,sgal.model_name)
            diag_loc = os.path.join(fig_loc,'diag')
            GenUtils.ensure_dir(diag_loc+'/')
            
            figname=os.path.join(diag_loc,'%s_%s_diag.png'%(ID,sgal.model_name))
            LFUtils.diagnostic_cmd(sgal,gal.trgb,figname=figname,inds=sgal.rel_ind)
            # should take fig_loc...
            LFUtils.plot_LFIR(gal,sgal,p_value,maglims)
            
            
            
        
        out.write('%s %s %.3f %i %i %i %i\n' %(sgal.ID,
                                               sgal.model,
                                               p_value,
                                               gal.irgb.size,
                                               gal.iagb.size,
                                               sgal.rel_rgb.size,
                                               sgal.rel_agb.size))
    
    out.close()    

    
    for model in models:
        compare_sims(IDs,model)
        if publish_plots:
            html_file = os.path.join(diag_loc,sgal.model_name+'_diag.html')
            gst.side_by_side(diag_loc,sgal.model_name,html_file)
            
            html_file = os.path.join(fig_loc,sgal.model_name+'.html')
            gst.one_col(fig_loc,sgal.model_name,html_file)
    
    return 




def compare_sims(IDs,model):
    sgals = []
    gals = []
    IDs = all_IDs()
    for ID in IDs:
        spread_file = os.path.join(model_src,'spread','spread_output_%s_model_%s'%(ID,model))
        #spread_file = os.path.join(model_src,'ast','ast_output_%s_model_%s'%(ID,model))
        gal = galaxy(ID,'ir')
        sgal = simgalaxy(spread_file,gal.filter1,gal.filter2)
        sgal.ID = ID
        sgal.model = model
        sgal.mix_modelname(model)
        maglims = (np.nan,np.nan)
        p_value = LFUtils.calc_LF(gal,sgal,maglims,normalize=False)
        sgals.append(sgal)
        gals.append(gal)

    # galaxies-ify
    SGals = galaxies(sgals)
    Gals = galaxies(gals)
    # select on the ones with measured z:
    galsz = Gals.finite_key('z')
    # sort by z.
    galsz_sort = sorted(galsz,key=lambda galaxy: galaxy.z)
    
    # use that to grab the sim gals
    sgalsz = []
    for g in galsz_sort:
        sgalsz += [sg for sg in SGals.galaxies if sg.target == g.target]
    
    zs = [g.z for g in galsz_sort]
    ax = plt.axes()
    for z,g,s in zip(zs,galsz_sort,sgalsz):
        ax.plot(z,float(s.rel_agb.size)/float(g.iagb.size),
                         'o',ms=5,color='black')
        bright_rgb = LFUtils.brighter(s.ast_mag2,g.trgb,inds = s.irgb)
        ax.plot(z,float((s.rel_agb.size-bright_rgb.size))/float(g.iagb.size),
                         'o',ms=5,color='red')        
        ax.annotate(g.target,xy=(z,float((s.rel_agb.size-bright_rgb.size))/float(g.iagb.size)+0.01),ha='center')
        ax.annotate(g.target,xy=(z,float(s.rel_agb.size)/float(g.iagb.size)+0.01),ha='center')
    ax.set_ylabel(r'$N_{AGB,model}/N_{AGB,data}$',fontsize=20)
    ax.set_xlabel(r'$Z$',fontsize=20)
    ax.set_title(r'$\rm{%s}$'%sgal.model_name.replace('_','\ '))
    GenUtils.ensure_dir(os.path.join(plt_dir,sgal.mix))
    fig_name = os.path.join(plt_dir,sgal.mix,sgal.model_name+'.png')
    plt.savefig(fig_name)
    plt.close()
    print 'wrote %s'%fig_name
    return sgal,gal
    
def write_spread_catalog(sgal,outfile=None):
    '''
    writes a slice of the trilegal output catalog.
    slice is on ast recovered stars that are randomly selected by normalization.
    i.e sgal.rec and sgal.rel_ind.
    returns string outfile name.
    '''
    if outfile==None:
        outfile = os.path.join(model_src,'spread',sgal.name.replace('ast','spread'))
    
    if not os.path.isfile(outfile):
        out = open(outfile,'w')
        # get the header
        sort_keys = sorted(sgal.data.key_dict.iteritems(), key=operator.itemgetter(1))
        col_keys =  [k[0] for k in sort_keys]
        
        # slice data on recovered inds and rel_inds. 
        data_slice = np.squeeze(sgal.data.data_array[[sgal.rec],:])[sgal.rel_ind]
        
        out.write('# {}\n'.format(' '.join(col_keys)))
        np.savetxt(out, data_slice, fmt='%-7.6f')
        print 'wrote %s'%outfile
    else:
        print '%s exists. Write spread catalog will not overwrite.'%outfile

    return outfile

if __name__=="__main__":
    #main(make_plots=True,publish_plots=True)
    
    # load all targets
    IDs = all_IDs()
    #IDs = [IDs[0]]
    
    # load all models
    models = ["cmd_input_CAF09_S_SCS.dat",
              "cmd_input_CAF09_S_SCSFG.dat",
              "cmd_input_CAF09_S_SCSFG_ETA2.dat"]
    
    for model in models:
        sgal,gal = compare_sims(IDs,model)
        fig_loc = os.path.join(plt_dir,sgal.mix,sgal.model_name)
        diag_loc = os.path.join(fig_loc,'diag')
        html_file = os.path.join(diag_loc,sgal.model_name+'_diag.html')
        gst.side_by_side(diag_loc,sgal.model_name,html_file)
        html_file = os.path.join(fig_loc,sgal.model_name+'.html')
        gst.one_col(fig_loc,sgal.model_name,html_file)