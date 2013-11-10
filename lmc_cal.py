import ResolvedStellarPops.graphics.GraphicsUtils as rspg
import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz
import numpy as np
import matplotlib.pyplot as plt
import os
import brewer2mpl
from TPAGBparams import research_path
import scipy.integrate
import galaxy_tests
from matplotlib.ticker import MultipleLocator
import sfh_tests
'''
actually for the LMC, you can simplify things a lot:

1- the distance and reddening are known, take them from Stefano's paper.
    you don't need to relocate the TRGB.

2- you don't need to renormalize the simulations to the RGB numbers, because
    that was already done by Stefano. Stefano provides the SFR in units of
    Msun/yr. Integrate this SFR(t) from 0 to 15 Gyr, and you have the total
    mass of stars ever formed in that region of the LMC. This total mass is
    given to TRILEGAL in the line:

1.0e6 10.0 # object_mass, object_dist: total mass inside field, distance

Provided that you use the same IMF and binary fraction as Stefano, you
shouldn't have any problem.

Doing so, the only data you actually need for the LMC is the 2MASS data
(which is 99% complete for the AGB stars), the SFR files Stefano gave you
already, and the values you find in Stefano's paper (distance, Av, etc).

Let me know in case this is not clear.
'''

class VaryVMCSFHs(sfh_tests.StarFormationHistories):
    def __init__(self, galaxy_input, filename='default', outfile_loc='default'):
        if filename == 'default':
            filename = research_path + 'TP-AGBcalib/LMC_Calib/SFR_LMC88.dat'
        self.prefix = os.path.split(filename)[1].split('.')[0]

        self.galaxy_input = galaxy_input
        if outfile_loc == 'default':
            self.outfile_loc = research_path + 'TP-AGBcalib/LMC_Calib/mc/'
        rsp.fileIO.ensure_dir(self.outfile_loc)

        self.data = self.parse_stefano_sfr()

        sfh_tests.StarFormationHistories.__init__(data=self.data)


    def prepare_outfiles(self):
        '''
        note, all attrs saved are file objects. Probably should have something
        to not keep adding to the same files on different runs. 
        '''
        # LF file
        lfname = os.path.join(self.outfile_loc, '%s_lf.dat' % self.prefix)
        self.lf_file = open(lfname, 'a')

        # mass_met_file
        mass_met_fname = os.path.join(self.outfile_loc,
                                      '%s_mass_met.dat' % self.prefix)

        self.mass_met_file = open(mass_met_fname, 'a')
        
        # CSLF file
        cslf_fname = os.path.join(self.outfile_loc,
                                  '%s_mass_met.dat' % self.prefix)
        self.cslf_file = open(cslf_fname, 'a')
        return
    
    def parse_stefano_sfr(self):
        orig_data = rsp.fileIO.readfile(self.filename, col_key_line=1)
    
        to = orig_data['Center_age_bin']
        half_lagebin = np.diff(to)
        tf = to[:-1] + half_lagebin
        tf = np.append(tf, 10.13)
    
        to = orig_data['Center_age_bin']
        to = np.insert(to, 0, 6.6)
        half_lagebin = np.diff(to)
        tf = to[:-1] + half_lagebin
        to = to[:-1]
    
        sfr = orig_data['med_SFR'] * 1e-3/ 2.
        sfr_errp = orig_data['max_SFR'] * 1e-3/ 2.
        sfr_errm = orig_data['min_SFR'] * 1e-3/ 2.
        
        feh = orig_data['med_FeH']
        feh_errp = orig_data['max_FeH']
        feh_errm = orig_data['min_FeH']
        
        mass = orig_data['med_MASS'] * 1e6
        mass_errp = orig_data['max_MASS'] * 1e6
        mass_errm = orig_data['min_MASS'] * 1e6
    
        dtype = [('lagei', '<f8'),
                 ('lagef', '<f8'),
                 ('sfr', '<f8'),
                 ('sfr_errp', '<f8'),
                 ('sfr_errm', '<f8'),
                 ('feh', '<f8'),
                 ('feh_errp', '<f8'),
                 ('feh_errm', '<f8'),
                 ('mass', '<f8'),
                 ('mass_errp', '<f8'),
                 ('mass_errm', '<f8')]
        
        # err... mh_disp is mass here and mh is feh...
        arr = np.vstack([to, tf, sfr, sfr_errp, sfr_errm, feh, feh_errp, feh_errm,
                         mass, mass_errp, mass_errm])
        data = arr.T.ravel().view(dtype)
        return data.view(np.recarray)
        
    
    def _make_trilegal_sfh(self, random_sfr=True, random_z=False,
                          zdisp=True, random_mass=True, outfile='default'):

        if outfile == 'default':
            outfile = os.path.join(self.base,
                                   self.name.replace('dat', '.sfr'))

        if random_z is False:
            feh = self.data.feh
        else:
            feh = self.random_draw_within_uncertainty('feh')

        z = convertz.convertz(feh=feh)[1]

        # this could be done better...
        # even dispersions
        zmin = convertz.convertz(feh=self.data.feh_errm)[1]
        zmax = convertz.convertz(feh=self.data.feh_errp)[1]
        zdisp = [np.mean([(zmax[i] - z[i]),
                          (z[i] - zmin[i])]) for i in range(len(z))]
        
        half_zbin = np.diff(z)/2.
        half_zbin = np.append(half_zbin, half_zbin[-1])
        z1 = z - half_zbin 
        z2 = z + half_zbin

        age1a = 10 ** (self.data.lagei)
        age1p = 1.0 * 10 ** (self.data.lagei + 0.0001)
        age2a = 1.0 * 10 ** self.data.lagef
        age2p = 1.0 * 10 ** (self.data.lagef + 0.0001)

        if random_sfr is False:
            sfr = self.data.sfr
        else:
            sfr = self.random_draw_within_uncertainty('sfr')

        if zdisp is True:
            zdisp = z * np.median(zdisp[np.nonzero(zdisp)])
            #zdisp = self.data.mh_disp
            fmt = '%.4e %.3e %.4f %.4f \n'
        else:
            zdisp = [''] * len(z)
            fmt = '%.4e %.3e %.4f %s\n'


        fmt = '%.4e %.3f %.4f %.4f\n'
        with open(outfile, 'w') as out:
            for i in range(len(sfr)):        
                out.write(fmt % (age1a[i], 0.0, z1[i], zdisp[i]))
                out.write(fmt % (age1p[i], sfr[i], z1[i], zdisp[i]))
                out.write(fmt % (age2a[i], sfr[i], z2[i], zdisp[i]))
                out.write(fmt % (age2p[i], 0.0, z2[i], zdisp[i]))
                out.write(fmt % (age1a[i], 0.0, z2[i], zdisp[i]))
                out.write(fmt % (age1p[i], sfr[i], z2[i], zdisp[i]))
                out.write(fmt % (age2a[i], sfr[i], z1[i], zdisp[i]))
                out.write(fmt % (age2p[i], 0.0, z1[i], zdisp[i]))
    
        mass = self.random_draw_within_uncertainty('mass')

        object_mass = scipy.integrate.simps(mass, self.data.lagei)
        return outfile, object_mass
    
    
    def prepare_trilegal(self, loidl=True, photsys='2mass', object_mass=None,
                         tri_inp_file='default', tri_sfr_file=None):
        tri_sfh_kw = tri_sfh_kw or {}

        if '2mass' in photsys:
            filter1 = 'Ks'
        elif 'ubv' in photsys:
            filter1 = 'K'

        if tri_inp_file == 'default':
            if loidl is True:
                tri_inp_file = tri_sfr_file.replace('.sfr','_loidl.inp')
            else:
                tri_inp_file = tri_sfr_file.replace('.sfr','.inp')
        
        gal_dict_inp = {'photsys': photsys,
                        'filter1': filter1,
                        'object_mass': object_mass,
                        'object_sfr_file': sfr_out_file}

        gal_dict = rsp.TrilegalUtils.galaxy_input_dict(**gal_dict_inp)
        gal_inp_pars = {'mag_limit_val': 20.6,
                        'binary_kind': 1,
                        'binary_frac': 0.3,
                        'object_sfr_mult_factorA': 1,
                        'object_dist': 47643.,
                        'object_av': 0.2, 
                        'object_avkind': 0}

        if loidl is True:
            gal_inp_pars['file_bcspec'] = ''
            gal_inp_pars['kind_mag'] = 2
        else:
            gal_inp_pars['kind_mag'] = 3
            if 'file_bcspec' in gal_inp_pars.keys():
                del gal_inp_pars['file_bcspec']

        gal_inp = rsp.fileIO.input_parameters(default_dict=gal_dict)
        gal_inp.add_params(gal_inp_pars)
        gal_inp.write_params(tri_inp_file,
                             rsp.TrilegalUtils.galaxy_input_fmt())
        return tri_inp_file
    
    def vary_the_SFH(self, cmd_input_file, prep_tri_kw=None, make_many_kw=None,
                     dry_run=False, diag_plots=False):
        
        tri_sfr_fmt = os.path.join(self.outfile_loc, 'tri_%s' % self.prefix)
        tri_sfr_fmt += '_003i.sfr'

        (sfr_out_file, object_mass) = zip(*[self._make_trilegal_sfh(outfile=tri_sfr_fmt % i, **make_many_kw)
                                            for i in range(nsfhs)])
        self.sfr_files = list(sfr_out_file)

        prep_tri_kw = dict({'loidl': True, photsys: '2mass'}.items() + prep_tri_kw.items())
        self.galaxy_inputs = [self.prepare_trilegal(tri_sfr_file=self.sfr_files[i], object_mass[i],
                                                    **prep_tri_kw)
                              for i in range(len(nsfhs))

        output = ADKLJSFL:KAJD
        # set up the trilegal run and then test!!!
        rsp.TrilegalUtils.run_trilegal(cmd_input_file, self.galaxy_input[i], output,
                                               loud=True)



def load_raw_vmc_data():
    #photom_file = research_path + 'TP-AGBcalib/LMC_Calib/photom.dat'
    #photom = rsp.fileIO.readfile(photom_file)
    results_file = research_path + 'TP-AGBcalib/LMC_Calib/photom_model.dat'
    dtype = [('RAJ2000', '<f8'), ('DEJ2000', '<f8'), ('recno', '<f8'),
             ('Seq', '<f8'), ('tau', '<f8'), ('logML', '<f8'), ('Lsun', '<f8'),
             ('Cl', '|S4'), ('Q', '|S1'), ('P', '|S1'), ('SED', '|S3'),
             ('recno1', '<f8'), ('Seq1', '<f8'), ('Umag', '<f8'),
             ('Bmag', '<f8'), ('Vmag', '<f8'), ('Imag', '<f8'),
             ('YmagV', '<f8'), ('JmagV', '<f8'), ('KsmagV', '<f8'),
             ('Jmag2', '<f8'), ('Hmag2', '<f8'), ('Ksmag2', '<f8'),
             ('[3.6]1', '<f8'), ('[3.6]2', '<f8'), ('[4.5]1', '<f8'),
             ('[4.5]2', '<f8'), ('[5.8]1', '<f8'), ('[5.8]2', '<f8'),
             ('[8.0]1', '<f8'), ('[8.0]2', '<f8'), ('[24]1', '<f8'),
             ('[24]2', '<f8'), ('chi2C', '<f8'), ('chi2O', '<f8')] 
    results = np.genfromtxt(results_file, dtype=dtype)
    return results

def read_vmc_table(filename):
    dtype = [('DEJ2000d', '<f8'),
             ('RAJ2000d', '<f8'),
             ('num', '<f8'),
             ('seq', '<f8'),
             ('RAJ2000', '|S12'),
             ('DEJ2000', '|S12'),
             ('tau', '<f8'),
             ('logMdot', '<f8'),
             ('Lum', '<f8'),
             ('Cl', '|S4'),
             ('Q', '|S1'),
             ('P', '|S1'),
             ('SED', '|S3'),
             ('Umag', '<f8'),
             ('Bmag', '<f8'),
             ('Vmag', '<f8'),
             ('Imag', '<f8'),
             ('YmagV', '<f8'),
             ('JmagV', '<f8'),
             ('KsmagV', '<f8'),
             ('Jmag2', '<f8'),
             ('Hmag2', '<f8'),
             ('Ksmag2', '<f8'),
             ('u3.6_1', '<f8'),
             ('u3.6_2', '<f8'),
             ('u4.5_1', '<f8'),
             ('u4.5_2', '<f8'),
             ('u5.8_1', '<f8'),
             ('u5.8_2', '<f8'),
             ('u8.0_1', '<f8'),
             ('u8.0_2', '<f8'),
             ('u24_1', '<f8'),
             ('u24_2', '<f8'),
             ('chi2C', '<f8'),
             ('chi2O', '<f8')]
    return np.genfromtxt(filename, dtype=dtype)


def read_lmc_cat(filename):
    with open(filename, 'r') as f:
        header = f.readline()
    col_keys = header.replace('#', '').strip().split()
    
    return np.genfromtxt(filename, names=col_keys)


def make_plot(output, extra='', photsys='2mass', tpagb_mass_bins=None):
    if tpagb_mass_bins is None:
        tpagb_mass_bins = np.arange(1., 6, 0.5)

    filter1 = 'J'
    if 'ubv' in photsys:
        filter2 = 'K'
    elif '2mass' in photsys:
        filter2 = 'Ks'
    sgal = rsp.Galaxies.simgalaxy(output, filter1=filter1, filter2=filter2)
    sgal.all_stages('TPAGB')
    sgal.mix_modelname(sgal.name)
    # http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/537/A105
    #vmcagbs = research_path + 'TP-AGBcalib/LMC_Calib/VMCAGBS_fmt.dat'
    #data = read_vmc_table(vmcagbs)

    vmcagbs = research_path + 'TP-AGBcalib/LMC_Calib/VMCAGBS_fmt2.dat'
    gal = rsp.Galaxies.galaxy(vmcagbs, filter1='Jmag2', filter2='Ksmag2', hla=False, 
                              angst=False)
    gal.filters = [filter1, filter2]
    #gal.maglims = [13., 12.]
    gal.target = sgal.mix

    smg = rsp.Galaxies.sim_and_gal(gal, sgal)
    outfile = output.replace('.dat','_%s.png' % extra)
    fig, axs, top_axs = smg.make_LF(filter1, filter2, color_hist=False,
                                    plot_tpagb=True , add_boxes=False,
                                    plot_LF_kw={'xlim': (0.5, 4)},
                                    figname=outfile)
    axs[1].cla()
    ax = sgal.color_by_arg(0,0,0, xdata=sgal.color, ydata=sgal.mag2,
                           coldata=sgal.data.get_col('m_ini'),
                           bins=tpagb_mass_bins, slice_inds=sgal.itpagb,
                           ax=axs[1], xlim=axs[0].get_xlim(),
                           ylim=axs[0].get_ylim(), fig=fig)

    tpagb_masses = sgal.data.get_col('m_ini')[sgal.itpagb]
    tpinds = np.digitize(tpagb_masses, tpagb_mass_bins)
    tpinds = tpinds[tpinds < len(tpagb_mass_bins)]
    tpagb_mass_inds = [sgal.itpagb[tpinds==i] for i in np.unique(tpinds)]
    # some masses are not recovered.
    tpagb_mass_bins = tpagb_mass_bins[np.unique(tpinds)]
    if 3 <= len(tpagb_mass_bins) <= 11:
        bmap = brewer2mpl.get_map('Paired', 'Qualitative', len(tpagb_mass_bins))
        cols = bmap.mpl_colors
    else:
        cols = rspg.discrete_colors(len(tpagb_mass_bins), colormap='RdYlGn')

    for i, tpagb_mass_ind in enumerate(tpagb_mass_inds):
        hist, _ = np.histogram(sgal.mag2[tpagb_mass_ind], bins=smg.bins)
        axs[2].semilogx(hist, smg.bins[1:], color=cols[i], ls='steps', lw=2)

    fig1, ax1 = plt.subplots()
    tpagb_imasses = sgal.data.get_col('m_ini')[sgal.itpagb]                
    tpagb_amasses = sgal.data.get_col('Mact')[sgal.itpagb]                
    dm = tpagb_imasses-tpagb_amasses
    ax1.plot(dm, sgal.mag2[sgal.itpagb], '.')
    
    outfile = outfile.replace('%s.png' % extra, '%s_by_mass.png' % extra)
    plt.savefig(outfile, dpi=300, bbox_to_inches='tight')
    print 'wrote %s' % outfile
    return sgal

def cslf(sgals, outfile=None):
    '''
    A simple CSLF plot.
    Takes a list of rsp.simgalaxy types.
    This takes the mag to be Mbol, and uses the first line of the trilegal
    catalogue to get dmod and Av.
    
    If an outfile is specified, will save fig.
    '''
    #sgals = [rsp.Galaxies.simgalaxy(tri_out, filter1='J', filter2='Ks',
    #                                photsys='2mass') for tri_out in tri_outs]

    [sgal.load_ic_mstar() for sgal in sgals]
    bmap = brewer2mpl.get_map('Set1', 'Qualitative', 3)
    cols = bmap.mpl_colors  
    bins = np.arange(-6.5, -2.5, 0.15)

    fig, ax = plt.subplots(figsize=(8, 8))
    for i, sgal in enumerate(sgals):
        model = sgal.name.split('_')[-1].split('.dat')[0]
        lab = galaxy_tests.translate_model_name(model) + ', $N=%i$' % (len(sgal.icstar))
        Mbol = sgal.data.get_col('mbol') - sgal.data.get_col('m-M0')[0] - sgal.data.get_col('Av')[0]
        plt.hist(Mbol[sgal.icstar], bins=bins, histtype='step', lw=5-i, color='white')
        plt.hist(Mbol[sgal.icstar], bins=bins, histtype='step', label=lab, lw=4-i, color=cols[i])
    ax.legend(frameon=False, prop={'size': 16}, loc=0)
    ax.set_xlim()[::-1]
    ax.set_xlabel('$M_{\\rm bol}$', fontsize=20)
    ax.set_ylabel('$N_c$', fontsize=20)
    plt.tick_params(labelsize=16)
    if outfile is not None:
        plt.savefig(outfile, dpi=300)
    return

def color_mdot_plot(sgals):
    nmodels = len(sgals)
    bmap = brewer2mpl.get_map('Set1', 'Qualitative', nmodels)
    cols = bmap.mpl_colors
    results = load_raw_vmc_data()
    data_color = results['Jmag2'] - results['Ksmag2']
    data_logML = results['logML']
    data_icstar = np.nonzero(results['Cl'] == 'C')
    data_imstar = np.nonzero(results['Cl'] == 'O')
    fig, (axs) = plt.subplots(nrows=nmodels, sharex=True, sharey=True, figsize=(8,8))
    for i, sgal in enumerate(sgals):
        model = sgal.name.split('_')[-1].split('.dat')[0]
        lab = galaxy_tests.translate_model_name(model)
        logml = sgal.data.get_col('logML')
        color = sgal.data.get_col('J') - sgal.data.get_col('Ks')
        axs[i].plot(color[sgal.icstar], logml[sgal.icstar], '.', label=lab, color=cols[i])
        axs[i].plot(color[sgal.imstar], logml[sgal.imstar], 'x', label=lab, color=cols[i])
        axs[i].plot(data_color[data_icstar], data_logML[data_icstar], '.', color='k')
        axs[i].plot(data_color[data_imstar], data_logML[data_imstar], 'x', color='k')
        axs[i].text(0.85, 0.1, lab, transform=axs[i].transAxes, fontsize=16)
    
    axs[1].set_ylabel('$\log \dot{M}\ ({\\rm M_\odot/yr})$', fontsize=20)
    axs[0].set_ylim(-11, -4)
    axs[0].yaxis.set_major_locator(MultipleLocator(2))
    axs[0].yaxis.set_minor_locator(MultipleLocator(.5))
    axs[0].xaxis.set_minor_locator(MultipleLocator(.2))
    axs[2].set_xlabel('$J-Ks$', fontsize=20)
    plt.tick_params(labelsize=16)
    fig.savefig('mdot_jks.png', dpi=300)
    return

def main():
    overwrite = True
    loidl = [True, False]
    photsys = '2mass'
    extra = ''
    for i in range(len(loidl)):
        if loidl is True:
            extra = '_loidl_%s' % photsys
        else:
            extra = '_%s' % photsys
        outputs = make_trilegal_sim(loidl=loidl[i],
                                    photsys=photsys,
                                    overwrite=overwrite,
                                    extra=extra)

        sgals = [make_plot(output, extra=extra) for output in outputs]
        outfile = os.path.join(os.path.split(output)[0], 'cslf%s.png' % extra)
        cslf(sgals, outfile=outfile)
        
if __name__ == "__main__":
    import pdb
    pdb.set_trace()
    main()