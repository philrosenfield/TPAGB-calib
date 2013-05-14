import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, NullFormatter
import os
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

def parse_stefano_sfr():    
    filename = '/Users/phil/research/TP-AGBcalib/LMC_Calib/SFR_LMC88.dat'
    outfile = '/Users/phil/research/TP-AGBcalib/LMC_Calib/tri_SFR_LMC88.dat'
    data = rsp.fileIO.readfile(filename, col_key_line=1)

    to = data['Center_age_bin']
    half_lagebin = np.diff(to)
    tf = to[:-1] + half_lagebin
    tf = np.append(tf, 10.13)

    sfr = data['med_SFR']

    z = convertz.convertz(feh=data['med_FeH'])[1]
    zmin = convertz.convertz(feh=data['min_FeH'])[1]
    zmax = convertz.convertz(feh=data['max_FeH'])[1]
    # even dispersions
    zdisp = [np.mean([(zmax[i] - z[i]),
                      (z[i] - zmin[i])]) for i in range(len(z))]

    half_zbin = np.diff(z)/2.
    half_zbin = np.append(half_zbin, half_zbin[-1])
    z1 = z - half_zbin 
    z2 = z + half_zbin

    age1a = 10 ** (to)
    age1p = 1.0 * 10 ** (to + 0.0001)
    age2a = 1.0 * 10 ** tf
    age2p = 1.0 * 10 ** (tf + 0.0001)

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
    
    object_mass = np.sum(data['med_MASS']) * 1e6
    return outfile, object_mass


def make_trilegal_sim(cmd_input=None, loidl=True, photsys='2mass', overwrite=True):

    if '2mass' in photsys:
        filter1 = 'Ks'
    elif 'ubv' in photsys:
        filter1 = 'K'

    cmd_inputs = ['/Users/phil/research/TP-AGBcalib/cmd_inputfiles/cmd_input_CAF09_S_MAR13.dat',
                  '/Users/phil/research/TP-AGBcalib/cmd_inputfiles/cmd_input_CAF09_S_APR13.dat']
    outputs = []
    for cmd_input in cmd_inputs:
        object_sfr_file, object_mass = parse_stefano_sfr()
        galaxy_input = object_sfr_file.replace('.dat','_galinp.dat')
        output = '%s_%s' % (object_sfr_file.replace('.dat',''),
                            cmd_input.split('cmd_input_')[1])
        
        gal_dict_inp = {'photsys': photsys,
                        'filter1': filter1,
                        'object_mass': object_mass,
                        'object_sfr_file': object_sfr_file}

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
            del gal_inp_pars['file_bcspec']

        gal_inp = rsp.fileIO.input_parameters(default_dict=gal_dict)
        gal_inp.add_params(gal_inp_pars)
        gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
        if os.path.isfile(output):
            if overwrite is True:
                rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input, output, loud=True)
            else:
                print output, 'found, not going to run trilegal and overwrite.'
        outputs.append(output)

    return outputs

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

def make_plot(output, extra='', photsys='2mass'):
    filter1 = 'J'
    if 'ubv' in photsys:
        filter2 = 'K'
    elif '2mass' in photsys:
        filter2 = 'Ks'
    sgal = rsp.Galaxies.simgalaxy(output, filter1=filter1, filter2=filter2)
    sgal.all_stages('TPAGB')
    sgal.mix_modelname(sgal.name)
    # http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/537/A105
    #vmcagbs = '/Users/phil/research/TP-AGBcalib/LMC_Calib/VMCAGBS_fmt.dat'
    #data = read_vmc_table(vmcagbs)

    vmcagbs = '/Users/phil/research/TP-AGBcalib/LMC_Calib/VMCAGBS_fmt2.dat'
    gal = rsp.Galaxies.galaxy(vmcagbs, filter1='Jmag2', filter2='Ksmag2', hla=False, 
                              angst=False)
    gal.filters = [filter1, filter2]
    #gal.maglims = [13., 12.]
    gal.target = sgal.mix

    smg = rsp.Galaxies.sim_and_gal(gal, sgal)
    outfile = output.replace('.dat','_mag%s.png' % extra)
    smg.make_LF(filter1, filter2, color_hist=True, plot_tpagb=True,
                add_boxes=False, plot_LF_kw={'xlim': (0.5, 4)}, figname=outfile)
    '''
    #fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
    #ax1.plot(lmc_cat['J_MAG'] - lmc_cat['Ks_MAG'], lmc_cat['Ks_MAG'],',', color='black')
    #ax2.plot(sgal.color, sgal.mag2, ',', color='black')
    #ax2.plot(sgal.color[sgal.itpagb], sgal.mag2[sgal.itpagb], '.', color='red')
    #ax2.set_ylim(ax2.get_ylim()[::-1])
    #ax1.set_xlim(ax2.get_xlim())
    #ax1.set_ylim(ax2.get_ylim())
    #ax1.set_xlabel('$J-Ks$')
    #ax1.set_ylabel('$Ks$')
    #ax2.set_title('$N_{TPAGB} = %i$' % len(sgal.itpagb))
    #[ax.tick_params(labelsize=16) for ax in [ax1, ax2]]
    
    
    # use the 2mass filters
    mag2 = data['Ksmag2']
    color = data['Jmag2'] - mag2


    fig, axs = plot_LF(color, mag2, sgal.color, sgal.mag2, sgal.filter1, sgal.filter2,
                       itpagb=sgal.itpagb, xlim=(-1, 2), ylim=(15,10),
                       model_plt_color='grey', outfile=outfile)

    # use the VMC filters
    extra += '_VMC'
    mag2 = data['KsmagV']
    color = data['JmagV'] - mag2

    outfile = output.replace('.dat','_mag%s.png' % extra)
    fig, axs = plot_LF(color, mag2, sgal.color, sgal.mag2, sgal.filter1, sgal.filter2,
                       itpagb=sgal.itpagb, xlim=(-1, 2), ylim=(15,10),
                       model_plt_color='grey', outfile=outfile)
    '''

def plot_LF(color, mag, scolor, smag, filt1, filt2,
            model_plt_color='red', data_plt_color='black', ylim=None,
            xlim=None, xlim2=None, model_title='Model', title=False,
            itpagb=None, res=0.1, outfile=None):

    def setup_lfplot(filt1, filt2, model_title='Model', target='LMC88',
                     lab_kw={}, data_plt_color='black',
                     model_plt_color='grey'):

        fig = plt.figure(figsize=(9, 9))
        # plot limits determined by hand
        bottom, height = 0.1, 0.8
        widths = [0.57, 0.22]
        lefts = [0.1, 0.73]

        axs = [plt.axes([lefts[i], bottom, widths[i], height])
               for i in range(len(lefts))]

        lab_kw = dict({'fontsize': 20}.items() + lab_kw.items())

        # titles
        #axs[0].set_title('$%s$' % target,
        #                 color=data_plt_color, **lab_kw)

        axs[0].set_title('$%s$' % model_title, color=model_plt_color, **lab_kw)

        axs[0].set_xlabel('$%s-%s$' % (filt1, filt2), **lab_kw)

        axs[0].set_ylabel('$%s$' % filt2, **lab_kw)
        #axs[1].set_xlabel(axs[0].get_xlabel(), **lab_kw)
        axs[1].set_xlabel('$\#$', **lab_kw)

        for ax in axs:
            ax.tick_params(labelsize=20)

        return fig, axs

    def fix_plot(axs, xlim=None, xlim2=None, ylim=None):
        # fix axes
        if xlim is not None:
            axs[0].set_xlim(xlim)
        if ylim is not None:
            axs[0].set_ylim(ylim)
        if xlim2 is not None:
            axs[1].set_xlim(xlim2)

        for ax in axs[:1]:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_minor_locator(MultipleLocator(0.2))

        for ax in axs:
            ax.set_ylim(axs[0].get_ylim())
            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))

        #axs[1].set_xlim(axs[0].get_xlim())

        # no formatters on mid and right plots
        [ax.yaxis.set_major_formatter(NullFormatter()) for ax in axs[1:]]


    fig, axs = setup_lfplot(filt1, filt2, model_title=model_title,
                            data_plt_color=data_plt_color,
                            model_plt_color=model_plt_color)

    axs[0].plot(color, mag, '.', color=data_plt_color)

    axs[0].plot(scolor, smag, '.', color=model_plt_color, mec=model_plt_color)
    if itpagb is not None:
        axs[0].plot(scolor[itpagb], smag[itpagb], '.', color='red')

    nbins = (mag.max() - mag.min()) / res

    gal_hist, bins = np.histogram(mag, nbins)
    sgal_hist, _ = np.histogram(smag, bins=bins)
    if itpagb is not None:
        sgal_hist2, _ = np.histogram(smag[itpagb], bins=bins)

    # plot histogram
    hist_kw = {'drawstyle': 'steps', 'color': data_plt_color, 'lw': 2}
    axs[1].semilogx(gal_hist, bins[1:], **hist_kw)


    if itpagb is not None:
        hist_kw['color'] = 'red'
        axs[1].semilogx(sgal_hist2, bins[1:], **hist_kw)
    else:
        hist_kw['color'] = model_plt_color
        axs[1].semilogx(sgal_hist, bins[1:], **hist_kw)

    fix_plot(axs, xlim=xlim, xlim2=xlim2, ylim=ylim)
    if outfile is not None:
        plt.savefig(outfile, dpi=300)
        print 'wrote %s' % outfile
    return fig, axs

def main():
    overwrite = True
    photsyss = ['2mass']
    for photsys in photsyss:
        outputs = make_trilegal_sim(loidl=True, photsys=photsys, overwrite=overwrite)
        extra = '_loidl_%s' % photsys
        [make_plot(output, extra=extra) for output in outputs]
        extra = '_%s' % photsys
        outputs = make_trilegal_sim(photsys=photsys, overwrite=overwrite)
        [make_plot(output) for output in outputs]

if __name__ == "__main__":
    import pdb
    pdb.set_trace()
    main()