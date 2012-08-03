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
        print 'reading tagged_file catalog'
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
    dict = {"SCL-DE1": ,
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
    return dict[ID]
'''
    
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
        
class simgalaxy(object):
    def __init__(self,trilegal_out,filter1,filter2):
        self.data = GenUtils.read_table(trilegal_out)
        self.diff1 = self.data.get_col('diff_'+filter1)
        self.diff2 = self.data.get_col('diff_'+filter2)
        self.recovered1 = np.nonzero(abs(self.diff1)<90.)[0]
        self.recovered2 = np.nonzero(abs(self.diff2)<90.)[0]
        self.rec = list(set(self.recovered1) & set(self.recovered2))
        self.mag1 = self.data.get_col(filter1)[self.rec]
        self.mag2 = self.data.get_col(filter2)[self.rec]
        self.stage = self.data.get_col('stage')[self.rec]
        self.color = self.mag1-self.mag2
        self.ast_mag1 = self.mag1+self.diff1[self.rec]
        self.ast_mag2 = self.mag2+self.diff2[self.rec]
        self.ast_color = self.ast_mag1-self.ast_mag2
        simgalaxy.load_ic_mstar(self)
        
    def stage_inds(self,stage_name):
        return np.nonzero(self.stage == get_stage_label(stage_name))[0]

    def load_ic_mstar(self):
        co = self.data.get_col('C/O')[self.rec]
        lage = self.data.get_col('logAge')[self.rec]
        mdot = self.data.get_col('logML')[self.rec]
        logl = self.data.get_col('logL')[self.rec]
        self.imstar, = np.nonzero((co<=1) & (logl >= 3.3) & (mdot<=-5) & (self.stage == get_stage_label('TPAGB')))
        self.icstar, = np.nonzero((co>=1) & (mdot<=-5) & (self.stage == get_stage_label('TPAGB')))
         
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

def main(ID,model):
    outfile = 'result_tab.dat'
    if os.path.isfile(outfile): 
        print outfile,'exists. appending.'
        out = open(outfile,'a')
    else:
        out = open(outfile,'w')
        out.write('# ID model p_value NRGB_data NAGB_data NRGB_model NAGB_model mass_model\n')

    # what's this for??
    band = 'ir'

    gal = galaxy(ID,band)

    # initializations
    object_mass = 5e6
    go = 0
    over_write = 0
    normalization = 1e9
    while normalization >.75:
        if go != 0: object_mass = object_mass*5.
        go +=1
        print 'Trying %s %s, Mass = %g, Attempt %i'%(ID,model,object_mass,go)
        trilegal_out = mk_sims.mk_sims(ID,model,object_mass=object_mass,over_write=over_write)
        #print 'trilegal_out has %i lines'%len(open(trilegal_out,'r').readlines())
        sgal = simgalaxy(trilegal_out,gal.filter1,gal.filter2)
        sgal.ID = ID
        sgal.model = model
        
        # how far below trgb?
        offset = 1.5
        maglims = (gal.trgb+offset,gal.trgb+0.5)
        p_value = LFUtils.calc_LF(gal,sgal,maglims)
        normalization = sgal.normalization
        over_write = 1
        print 'normalization',normalization

    print 'necessary object_mass = %g'%object_mass
    out.write('%s %s %.3f %i %i %i %i %e\n' %(ID,model,p_value,gal.irgb.size,gal.iagb.size,sgal.rel_rgb.size,sgal.rel_agb.size,object_mass))
    LFUtils.diagnostic_cmd(sgal,gal.trgb,figname=os.path.join(plt_dir,'%s_%s_diag.png'%(ID,model)))
    LFUtils.plot_LFIR(gal,sgal,p_value,maglims)
    out.close()


if __name__=="__main__":
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
        "DDO71"
        ]
    
    models = ["cmd_input_CAF09_S_SCS.dat",
              "cmd_input_CAF09_S_SCSFG.dat",
              "cmd_input_CAF09_S_SCSFG_ETA2.dat"]
    

    
        

    #for ID in IDs: mk_sims.mk_sims(ID,models)
    #sys.exit()
    #run_all(IDs,models)
    for ID,model in itertools.product(IDs,models):
        main(ID,model)