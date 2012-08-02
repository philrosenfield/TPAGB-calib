import sys,os
#from PadovaTracksUtils import spread_angst
from subprocess import *

import numpy as np
from TPAGBparams import *

photosys="wfc3snap"
mag_num=9

def get_sfrFILE(ID):
    sfr_dir = os.path.join(data_src,'sfh')
    l=os.listdir(sfr_dir)
    IDs=[f.split("_")[2] for f in l if not f.startswith('.')]
    ind=IDs.index(ID)
    return os.path.join(sfr_dir,l[ind])

def get_fakFILE(ID):
    fak_dir =os.path.join(data_src,'fakes')
    l=os.listdir(fak_dir)
    IDs=[f.split("_")[1] for f in l if not f.startswith('.')]
    if ID.startswith("UGC-"):
        tmp=ID.split("-")
        ID="%s%s"%(tmp[0],int(tmp[1]))
        if len(tmp)>2:
            ID+="-"+"-".join(tmp[2:])

    if ID=="NGC0404-DEEP": ID="NGC404"
    ind=IDs.index(ID)
    return os.path.join(fak_dir,l[ind])


def mk_sims(ID,
            models,
            mfactor=1.,
            outdir=model_src,
            object_mass=1e6,
            over_write=False):
    #print ID
    # read data from table and select parameters
    # for galaxy ID
    trilegal_outs = []
    table=open(os.path.join(table_src,'table.dat'))
    ok=0
    for lines in table.readlines():
        line=lines.strip().split()
        if line[0]==ID:
            (dmod,
             object_av,
             camera,
             filtname,
             maglim,
             )=line[1:]
            ok=1
            #if camera=="ACS":
            #    mag_num=1
            #    mag_lim=31
            #elif camera=="WFPC2":
            #    mag_num=1
            #    mag_lim=31
            #else:
            #    print "WRONG CAMERA NAME!"
            #    sys.exit(2)
                
    table.close()
    if not ok:
        print "WRONG ID"
        sys.exit(2)
 
    object_dist = 10**((5+float(dmod))/5)

    # define filenames
    #sfr_file="sfr_%s.dat"%ID
    sfr_file=get_sfrFILE(ID)

    if (os.path.isfile(sfr_file)):
        object_sfr=os.path.abspath(sfr_file)
    else:
        print "no such file %s"%sfr_file
        sys.exit(2)
    
    def_file    = os.path.join(model_src,'pars','default.pars') # file with default params

    par_dir="%s/pars"%outdir
    par_filename    ="run_trilegal_input_%s.pars"%ID # run_trilegal par file
    par_file="%s/%s"%(par_dir,par_filename)

    inp_dir="%s/input"%outdir
    inp_filename    ="input_%s.dat"%ID    # trilegal input file
    inp_file="%s/%s"%(inp_dir,inp_filename)

    if not os.path.isdir(par_dir):
        os.makedirs(par_dir)
    if not os.path.isdir(inp_dir):
        os.makedirs(inp_dir)
    
    fak_file=get_fakFILE(ID)

    mag_lim=np.loadtxt(fak_file,usecols=(1,)).max()+2
    if not os.path.isfile(fak_file):
        print "no such file %s"%fak_file
        sys.exit(2)

    ifile=open(def_file)
    pfile=open(par_file,'w')

    # write params file for run_trilegal
    fmt="%-18s %s\n"
    lines=ifile.read()
    pfile.write(lines)

    pfile.write(fmt%('photosys    ',photosys    )) 
    pfile.write(fmt%('mag_num     ',mag_num    )) 
    pfile.write(fmt%('mag_lim     ',mag_lim    )) 
    pfile.write(fmt%('object_mass ',object_mass)) 
    pfile.write(fmt%('object_dist ',object_dist)) 
    pfile.write(fmt%('object_av   ',object_av  )) 
    pfile.write(fmt%('object_sfr  ',object_sfr )) 
    
    pfile.close()
    ifile.close()
    
    # run trilegal
    #for model in models:
    # EDIT - made it only one model at a time...
    model = models
    #
    out_filename="output_%s_model_%s"%(ID,model)
    out_dir="%s/output"%outdir
    ast_filename="ast_%s"%out_filename
    ast_dir="%s/ast"%outdir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(ast_dir):
        os.makedirs(ast_dir)
    out_file    ="%s/%s"%(out_dir,out_filename)   # trilegal output file
    ast_file    ="%s/%s"%(ast_dir,ast_filename)   # trilegal output file + completeness + errors
    
    #cmd="$HOME/SOFTWARE/WXTRILEGAL/run_trilegal.py "
    cmd="/Users/phil/research/PyTRILEGAL/run_trilegal.py "
    #cmd="run_trilegal.py "
    cmd+="-e code/main "
    cmd+="-a "
    cmd+="-l "
    cmd+="-i %s "%os.path.abspath(inp_file)
    cmd+="-o %s "%os.path.abspath(out_file)
    if not model.startswith('cmd'): model = 'cmd_input_'+model
    cmd+="-f ../cmd_inputfiles/%s "%model
    cmd+=par_file
    
    if not os.path.isfile(out_file) or over_write==True:
        print 'running TRILEGAL:',model,ID
        print cmd
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,close_fds=True)
        stdout,stderr = (p.stdout,p.stderr)
        p.wait()
        
        EOF = os.path.join(os.environ['PYTHONPATH'],'EOF')
        cmd= '%s << %s\n'%(os.path.join(model_src,'spread_angst'),EOF)
        cmd+=os.path.abspath(fak_file)+"\n"
        cmd+=os.path.abspath(out_file)+"\n"
        cmd+=os.path.abspath(ast_file)+"\n"
        cmd+="%s\n"%EOF
        print "  ... completeness using %s"%fak_file
        print "  %s -> %s"%(out_file,ast_file)
        print 'Running spread_angst...'
        #p = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE)
        #sts = os.waitpid(p.pid, 0)[1]
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,close_fds=True)
        stdout,stderr = (p.stdout,p.stderr)
        p.wait()
    
        os.system("wc -l %s %s|head -2"%(out_file,ast_file))
    
    return os.path.abspath(ast_file)
    