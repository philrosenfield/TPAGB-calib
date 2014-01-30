import sfh_tests_multi_proc
import numpy as np
import ResolvedStellarPops as rsp
import galaxy_tests
import os
import tables
import matplotlib.pylab as plt
from TPAGBparams import snap_src
import logging
logger = logging.getLogger()
angst_data = rsp.angst_tables.AngstTables()

class StatisticalComparisons(object):
    def __init__(self, cmd_input_file, target, outfile_loc='default',
                 extra_str='', mc=True):
        self.target = target
        self.outfile_loc, self.fnames, self.agb_mod = \
            sfh_tests_multi_proc.setup_files(cmd_input_file, target,
                                             outfile_loc=outfile_loc,
                                             extra_str=extra_str, mc=mc)
        self.files = sfh_tests_multi_proc.FileIO()
        self.files.mc = mc

    def poission_chi2(self, hist_it_up=False, table_file='default'):
        self.files.ags = sfh_tests_multi_proc.load_default_ancient_galaxies(table_file=table_file)
        self.files.load_data_for_normalization(target=self.target, ags=self.files.ags)

        opt_gal, ir_gal = self.files.load_galaxies(hist_it_up=hist_it_up)

        # cut LF at 90% completeness
        obins, = np.nonzero(opt_gal.bins <= self.files.opt_offset)
        ibins, = np.nonzero(ir_gal.bins <= self.files.ir_offset)

        agb_obins, = np.nonzero(opt_gal.bins <= self.files.opt_trgb - \
                                                self.files.opt_trgb_err * \
                                                self.files.ags.factor[0])
        agb_ibins, = np.nonzero(ir_gal.bins <= self.files.ir_trgb - \
                                                self.files.ir_trgb_err * \
                                                self.files.ags.factor[1])
        opt_model_hists, opt_models_binss = self.files.load_lf_file(self.fnames[0])
        ir_model_hists, ir_models_binss = self.files.load_lf_file(self.fnames[1])
        opt_chi2 = np.array([])
        ir_chi2 = np.array([])

        opt_chi2_agb = np.array([])
        ir_chi2_agb = np.array([])
        nhists = np.min([len(ir_model_hists), len(opt_model_hists)])
        for i in range(nhists):
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(opt_gal.hist[obins[1:]],
                                                           opt_model_hists[i][obins[1:]])
            opt_chi2 = np.append(opt_chi2, chi2)
            #print 'opt', chi2, np.mean(sig)
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(ir_gal.hist[ibins[1:]],
                                                           ir_model_hists[i][ibins[1:]])
            ir_chi2 = np.append(ir_chi2, chi2)
            #print 'ir ', chi2, np.mean(sig)

            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(opt_gal.hist[agb_obins[1:]],
                                                           opt_model_hists[i][agb_obins[1:]])
            opt_chi2_agb = np.append(opt_chi2_agb, chi2)
            #print 'opt', chi2, np.mean(sig)
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(ir_gal.hist[agb_ibins[1:]],
                                                           ir_model_hists[i][agb_ibins[1:]])
            ir_chi2_agb = np.append(ir_chi2_agb, chi2)
            #print 'ir ', chi2, np.mean(sig)
        return opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb


def result2dict(result_files):
    res_dict = {}

    if 'narratio' in result_files[0]:
        search = 'ratio'
        chi2 = False
    elif 'chi2' in result_files[0]:
        search = 'chi2'
        chi2 = True
    else:
        print 'either narratio file or chi2 file'
        return {}

    for result_file in result_files:
        data = rsp.fileIO.readfile(result_file, string_column=0)
        target = os.path.split(result_file)[1].split('_')[3]
        agb_mod = os.path.split(result_file)[1].split(target)[0][:-1]
        fields = [d for d in data.dtype.names if search in d]
        key = '%s_%s' % (target, agb_mod)
        for f in fields:
            res_dict['%s_%s_mean' % (key, f)] = np.mean(data[f])
            if chi2 is True:
                res_dict['%s_%s_std' % (key, f)] = np.std(data[f])

    # to get narratio means
    #opt_total = [v for k,v in narr_dict.items() if 'opt_ar_ratio_mean' in k]
    #ir_total = [v for k,v in narr_dict.items() if 'ir_ar_ratio_mean' in k]
    #narr_dict['%s_opt_ratio_mean' % agb_mod] = np.mean(opt_total)
    #narr_dict['%s_ir__ratio_mean' % agb_mod] = np.mean(ir_total)
    #print '%s opt $%.2f\pm%.2f$' % (self.agb_mod, np.mean(opt_total), np.std(opt_total))
    #print '%s ir $%.2f\pm%.2f$' % (self.agb_mod, np.mean(ir_total), np.std(ir_total))
    # to get chi2 means
    #agb_mods = np.unique(agb_mods)
    #for extra in extras:
    #    for band in bands:
    #        band += extra
    #        for agb_mod in agb_mods:
    #            total = [v for k,v in chi_dict.items() if '%s_%s_mean' % (agb_mod, band) in k]
    #            chi_dict['%s_%s_mean' % (agb_mod, band)] = np.mean(total)
    #            chi_dict['%s_%s_std' % (agb_mod, band)] = np.std(total)

    return res_dict


def write_chi2_table(targets, cmd_input_files, table_file='default',
                     outfile_loc='default'):
    chi2_files = []
    for target in targets:
        for cmd_input_file in cmd_input_files:
            st = StatisticalComparisons(cmd_input_file, target,
                                        outfile_loc=outfile_loc)
            opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb = \
                st.poission_chi2(table_file=table_file)

            chi2_file = os.path.join(outfile_loc,
                                     '%s_%s_chi2.dat' % (st.agb_mod, target))
            cfmt = '%i %.3f %.3f %.3f %.3f \n'
            with open(chi2_file, 'w') as c2:
                c2.write('# sfr opt_chi2 ir_chi2 opt_agb_chi2 ir_agb_chi2\n')
                for i in range(len(opt_chi2)):
                    c2.write(cfmt % (i, opt_chi2[i], ir_chi2[i],
                                     opt_chi2_agb[i], ir_chi2_agb[i]))
            chi2_files.append(chi2_file)
    return chi2_files


def get_data(table_file='default'):
    if table_file == 'default':
        table_file = snap_src + '/tables/ancients_0.1_0.2_galaxies.dat'
    ags = sfh_tests_multi_proc.AncientGalaxies()
    ags.read_trgb_table(table_file)
    data_dict = {}
    for i, target in enumerate(ags.data.target):
        for band in ['opt', 'ir']:
            ratio = ags.data[i]['n%s_agb' % band] / \
                    ags.data[i]['n%s_rgb' % band]
            data_dict['%s_%s' % (target, band)] = ratio
            data_dict['%s_%s_err' % (target, band)] = \
                galaxy_tests.count_uncert_ratio(ags.data[i]['n%s_agb' % band],
                                                ags.data[i]['n%s_rgb' % band])
    return data_dict


def narratio_table(narratio_files, table_file='default'):
    '''write the latex table'''
    data_dict = get_data(table_file=table_file)
    nar_dict = result2dict(narratio_files)
    targets = list(np.unique([k.split('_')[0] for k in nar_dict.keys()]))
    agb_mods = list(np.unique(['_'.join(k.split('_')[1:4])
                               for k in nar_dict.keys()]))
    # table columns: target, data ratio, (ratio, frac diff) per agb mod
    ir_table = np.empty(shape=(len(targets), len(agb_mods)*2 + 2), dtype='|S20')
    opt_table = np.empty(shape=(len(targets), len(agb_mods)*2 + 2), dtype='|S20')
    fmt = r'$%.3f\pm%.3f$ & '
    fmt2 = r'$%.3f\pm%.3f$ \\'
    for key, val in nar_dict.items():
        # go through half the dict (see below)
        if 'err' in key:
            continue

        # choose which table
        if 'ir' in key:
            table = ir_table
            band = 'ir'

        if 'opt' in key:
            table = opt_table
            band = 'opt'

        # choose the correct row and column placement in the table
        target = key.split('_')[0]
        agb_mod = '_'.join(key.split('_')[1:4])
        row = targets.index(target)
        # column 0 is target, column 1 is data columns 2, 4, 6 have ratios
        column = (agb_mods.index(agb_mod) + 1) * 2

        # target
        table[row, 0] = '%s &' % target.upper()

        # data
        dnarr = data_dict['%s_%s' % (target.lower(), band)]
        derr = data_dict['%s_%s_err' % (target.lower(), band)]
        dstr = fmt % (dnarr, derr)
        table[row, 1] = dstr

        # model
        mnarr = val
        # grab the error from the dict
        err_key = key.replace('mean', 'err_mean')
        mnerr = nar_dict[err_key]
        mstr = fmt % (mnarr, mnerr)
        table[row, column] = mstr

        # frac difference
        pct_diff = (mnarr - dnarr) / dnarr
        pct_diff_err = np.abs(pct_diff * (mnerr/mnarr + derr/dnarr))
        f = fmt
        # if final column, put \\ not &
        if column == len(agb_mods) * 2:
            f = fmt2
        pdstr = f % (pct_diff, pct_diff_err)
        table[row, column + 1] = pdstr

    # write the file
    ratio_str = '$\\frac{N_{\\rm TP-AGB}}{N_{\\rm RGB}}$'
    header = 'Target & %s Data & ' % ratio_str
    for i in range(len(agb_mods)):
        header += '%s %s & Frac. Difference & ' % (ratio_str,
                                                   agb_mods[i].split('_')[-1])
    header = header[:-2] + r'\\ \n \hline \n'

    outfile_dir = os.path.split(narratio_files[0])[0]
    for band, table in zip(['opt', 'ir'], [opt_table, ir_table]):
        outfile = os.path.join(outfile_dir, '%s_narratio_table.tex' % band)
        with open(outfile, 'w') as out:
            out.write(header)
            np.savetxt(out, table, fmt='%s')


def data_table(targets):
    # target, av, dist, [opt: frac comp, trgb, nrgb, nagb ratio] [ir: ..]
    table_file = snap_src + '/tables/ancients_0.1_0.2_galaxies.dat'
    ags = sfh_tests_multi_proc.AncientGalaxies()
    ags.read_trgb_table(table_file)
    data_dict = get_data(table_file=table_file)
    comp_data = tables.read_completeness_table()
    row = ''
    row2 = ''
    opt_agb_tot = np.sum(ags.data.nopt_agb)
    opt_rgb_tot = np.sum(ags.data.nopt_rgb)
    ir_agb_tot = np.sum(ags.data.nir_agb)
    ir_rgb_tot = np.sum(ags.data.nir_rgb)

    totfmt = 'Total & %i & %i & $%.3f\\pm%.3f$ &  %i & %i & $%.3f\\pm%.3f$ \\\\'
    total = totfmt % (opt_agb_tot, opt_rgb_tot, opt_agb_tot/opt_rgb_tot,
                      galaxy_tests.count_uncert_ratio(opt_agb_tot, opt_rgb_tot),
                      ir_agb_tot, ir_rgb_tot, ir_agb_tot/ir_rgb_tot,
                      galaxy_tests.count_uncert_ratio(ir_agb_tot, ir_rgb_tot))

    for target in targets:
        if target == 'NGC2976-DEEP':
            extra_key='F606W,F814W'
        else:
            extra_key=None
        (Av, dmod) = [angst_data.get_item(target, i, extra_key=extra_key) for i in ['Av', 'dmod']]
        comp_row = rsp.fileIO.get_row(comp_data, 'target', target)
        nstars_row = rsp.fileIO.get_row(ags.data, 'target', target)
        sub_dict  = {}
        opt_err_pct = []
        ir_err_pct = []
        for (k,v) in data_dict.items():
            if '404' in target:
                target = target.lower().replace('-deep', '')
            if target in k or target.lower() in k:
                sub_dict[k] = v
        opt_err_pct.append(sub_dict['%s_opt_err' % target.lower()] /
                           sub_dict['%s_opt' % target.lower()])
        ir_err_pct.append(sub_dict['%s_ir_err' % target.lower()] /
                          sub_dict['%s_ir' % target.lower()])

        row += '%s & %.2f & %.2f & ' % (target, Av, dmod)
        row += '%(opt_filter2).2f & ' % comp_row
        row += '%(opt_trgb).2f & ' % nstars_row
        row += '%(ir_filter2).2f & ' % comp_row
        row += '%(ir_trgb).2f \\\\ \n' % nstars_row

        row2 += '%s & ' % target
        row2 += '%(nopt_agb)i & %(nopt_rgb)i & ' % nstars_row
        row2 += '$%.3f\\pm%.3f$ & ' % (sub_dict['%s_opt' % target.lower()], sub_dict['%s_opt_err' % target.lower()])
        row2 += '%(nir_agb)i & %(nir_rgb)i & ' % nstars_row
        row2 += '$%.3f\\pm%.3f$ \\\\ \n' % (sub_dict['%s_ir' % target.lower()], sub_dict['%s_ir_err' % target.lower()])

    print row
    print
    print row2
    print total

    print np.max(opt_err_pct), np.max(ir_err_pct)


def chi2plot(model_dict, outfile_loc=None):
    targets = list(np.unique([k.split('_')[0] for k in model_dict.keys()]))
    agb_mods = list(np.unique(['_'.join(k.split('_')[1:4])
                               for k in model_dict.keys()]))

    opt_filter = '$F814W$'
    ir_filter = '$F160W$'
    agb_sym = 's'
    tot_sym = 'o'
    opt_offset = -0.5
    ir_offset = 0.5
    cols = ['black', 'navy', 'darkred', 'darkgreen', 'purple', 'orange']

    fig, axs = plt.subplots(ncols=3)
    ymaxs = np.zeros(len(axs)) - 1.

    for key, val in model_dict.items():
        # color = target
        # symbol = agb or total
        if 'std' in key:
            continue
        target = key.split('_')[0]
        agb_mod = '_'.join(key.split('_')[1:4])
        ax_num = agb_mods.index(agb_mod)
        ax = axs[ax_num]
        col = cols[targets.index(target)]
        sym = tot_sym
        if '_agb' in key:
            sym = agb_sym
        offset = opt_offset
        if 'ir' in key:
            offset = ir_offset
        ax.plot(offset, val, sym, color=col, ms=20, alpha=0.6)
        ymaxs[ax_num] = np.max([val, ymaxs[ax_num]])

    axs[0].annotate('$\chi^2$', (.05, .5), fontsize=20, rotation='vertical',
                    xycoords='figure fraction')

    for i, ax in enumerate(axs):
        ax.set_ylim(0, ymaxs[i])
        ann_kw = {'ha': 'center', 'va': 'bottom', 'fontsize': 16}
        ax.annotate(opt_filter, (opt_offset, ymaxs[i]), **ann_kw)
        ax.annotate(ir_filter, (ir_offset, ymaxs[i]), **ann_kw)
        ax.set_xlabel('$%s$' % agb_mods[i].split('_')[-1].upper(),
                      fontsize=20)
        ax.set_xlim(-1, 1)
        ax.tick_params(labelsize=20, bottom='off', top='off', right='off',
                       labelbottom='off')

    # fake the legend
    ax.plot(-99, -99, agb_sym, ms=12, alpha=0.6, color='black',
            label='$AGB\ Only$')
    [ax.plot(-99, 99, tot_sym, ms=12, alpha=0.6, color=cols[j],
             label='$%s$' % targets[j].upper().replace('-DEEP', ''))
     for j in range(len(targets))]
    axs[-1].legend(loc='upper right', numpoints=1, bbox_to_anchor=(1.3, .95),
                   frameon=False)

    if outfile_loc is None:
        outfile_loc = os.getcwd()
    outfile = os.path.join(outfile_loc, 'chi2_plot.png')
    plt.savefig(outfile, dpi=150)
    return axs


'''
I like the plot more.
def chi2_table(targets, cmd_input_files, table_file='default', outfile_loc='default'):
    targets = galaxy_tests.load_targets(targets)
    fmt = '$%.3f\\pm%.3f$ & '
    bands = ['opt', 'ir']
    extras = ['', '_agb']
    model_dicts = []
    agb_mod = []
    for i in range(len(cmd_input_files)):
        write_chi2_table(targets, cmd_input_files, table_file='default', outfile_loc=outfile_loc)

        model_dict = result2dict(targets)

        header = 'Target & '
        for extra in extras:
            for band in bands:
                band += extra
                agb_mod_short = agb_mod.split('_')[-1]
                header += '\\chi^2 %s %s & ' % (agb_mod_short, band)
        header = header[:-2] + '\\\\'
        print header
        for target in targets:
            line = '%s & ' % target.upper()
            for extra in extras:
                for band in bands:
                    band += extra
                    line += fmt % (model_dict['%s_%s_mean' % (target, band)],
                                   model_dict['%s_%s_std' % (target, band)])
            print line
        model_dicts.append(model_dict)

        line = 'Mean $\\chi^2$ & '
        for extra in extras:
            for band in bands:
                band += extra
                line += fmt % (model_dict['%s_%s_mean' % (agb_mod, band)],
                               model_dict['%s_%s_std' % (agb_mod, band)])
        print line
    return model_dicts

'''
