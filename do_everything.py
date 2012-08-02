from LFUtils import *
from mk_sims import *
from multiprocessing import Pool
import itertools
import time

def run_all(IDs,models):
    pool=Pool()
    
    res=[]
    for ID,model in itertools.product(IDs,models):
        
        res.append(pool.apply_async(main,(ID,model),))

    for r in res:
        r.get()
    
    return

def main(ID,model):
    fits_src = os.path.join(data_src,'galaxies')
    outfile = 'result_tab.dat'
    if os.path.isfile(outfile): 
        print outfile,'exists. appending.'
    # what's this for??
    band = 'ir'
    out = open(outfile,'a')
    out.write('# ID model p_value NRGB_data NAGB_data NRGB_model NAGB_model mass_model\n')
    # add fluxs -- check with Jason's paper.
    
    if band == 'opt':
        fits = get_afile(fits_src,'*'+'*'.join((ID,'trim','.fits')))[0]
        trgb = get_tab5_trgb_Av_dmod(ID)[0]
    elif band == 'ir':
        # read data
        if ID == "NGC0404-DEEP": ID = "NGC404"
        fits, = get_afile(fits_src,'*'+'*'.join((ID,'IR','.fits')))
        trgb = get_trgb_ir_nAGB(ID)[0]
    else: 
        print 'choose opt or ir'
        sys.exit()
    
    f = readbintab(fits)
    filt1 = f['names'][-1].split('.')[0].split('_')[-2]
    filt2 = f['names'][-1].split('.')[0].split('_')[-1]
    print 'fitstable =',f['names'][-1]
    mag1 = f['Mag1']
    mag2 = f['Mag2']
    color = mag1-mag2

    if ID == "NGC0404-DEEP": ID = "NGC404"
    object_mass = 1e7
    go = 1
    normalization = 1e9
    #p_value,norm,nB_AGB,ps_nB_AGB = plot_LFIR(ID,model)
    while normalization >1.:
        print 'Trying',ID,model,'Mass =',object_mass,'try number...........',go
        trilegal_out = mk_sims(ID,model,object_mass=object_mass)
        
        synthcmd = read_table(trilegal_out)
        s_mag2 = synthcmd.get_col(filt2) + synthcmd.get_col('diff_'+filt2.strip())
        s_mag2 = s_mag2[np.nonzero(abs(synthcmd.get_col('diff_'+filt2.strip())) < 90.)[0]]

        Norm = trgb + 1.5
        # ind,nB_AGB,nNorm,ps_nNorm,ps_nB_AGBm,hist,bins,s_hist_normed,p_value,normalization
        calc_LF_out = calc_LF(mag2,s_mag2,Norm,trgb)
        ind,nB_AGB,nNorm,ps_nNorm,ps_nB_AGBm,hist,bins,s_hist_normed,p_value,normalization = calc_LF_out
        object_mass = object_mass*5.
        go +=1
        
    out.write('%s %s %.3f %i %i %i %i %e\n' %(ID,model,p_value,nNorm,nB_AGB,ps_nNorm,ps_nB_AGBm,object_mass))
    plot_LFIR(ID,model,trgb,Norm,filt1,filt2,color,mag2,synthcmd,*calc_LF_out)
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
        #"NGC0404-DEEP", # This is the optical name...
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

    models=["cmd_input_CAF09_S_SCS.dat"]
    '''
        #"cmd_input_ma08.dat",
        "cmd_input_gi10_rev.dat"]#,
        #"cmd_input_agb_wv.dat"        
        #]'''
    
    #for ID in IDs: mk_sims(ID,models)
    #sys.exit()
    for ID,model in itertools.product(IDs,models):
        main(ID,model)