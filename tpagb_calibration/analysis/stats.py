import logging
import os

import numpy as np
import matplotlib.pyplt as plt
import ResolvedStellarPops as rsp

from TPAGBparams import snap_src

logger = logging.getLogger()

angst_data = rsp.angst_tables.angst_table.AngstTables()

def narratio_table(self):
    narratio_files = rsp.fileIO.get_files(self.outfile_dir, '*narratio*dat')
    stats.narratio_table(narratio_files)
    return

def chi2_stats(targets, cmd_inputs, outfile_dir='default', extra_str=''):
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                            outfile_loc=outfile_dir,
                                            extra_str=extra_str)
    chi2_dicts = stats.result2dict(chi2_files)
    stats.chi2plot(chi2_dicts, outfile_loc=outfile_dir)
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                            outfile_loc=outfile_dir,
                                            extra_str=extra_str,
                                            just_gauss=True)
    return
def contamination_files(filenames):
    opt_eagb_contam = np.array([])
    opt_rheb_contam = np.array([])
    ir_eagb_contam = np.array([])
    ir_rheb_contam = np.array([])
    opt_ms_contam = np.array([])
    opt_bheb_contam = np.array([])
    ir_ms_contam = np.array([])
    ir_bheb_contam = np.array([])

    if type(filenames) == str:
        filenames = list(filenames)
    for filename in filenames:
        with open(filename, 'r') as fhandle:
            lines = fhandle.readlines()
        # rc contamination 12::13
        rgb_opt = [l for l  in lines if l.startswith('rgb opt')]
        rgb_data = zip(*[t.strip().split()[2:] for t in rgb_opt])
        rgb_data = np.array(rgb_data, dtype=float)
        eagb_in_rgb = rgb_data[5]/rgb_data[7]
        rheb_in_rgb = rgb_data[4]/rgb_data[7]
        opt_eagb_contam = np.append(opt_eagb_contam, np.max(eagb_in_rgb))
        opt_rheb_contam = np.append(opt_rheb_contam, np.max(rheb_in_rgb))
        #print filename, 'opt', np.max(eagb_in_rgb), np.max(rheb_in_rgb)

        opt =  [l for l  in lines if l.startswith('rgb opt') or l.startswith('agb opt')]
        data = zip(*[t.strip().split()[2:] for t in opt])
        data = np.array(data, dtype=float)
        ms_in_opt = data[0]/data[7]
        bheb_in_opt = data[3]/data[7]
        opt_bheb_contam = np.append(opt_bheb_contam, np.max(bheb_in_opt))
        opt_ms_contam = np.append(opt_ms_contam, np.max(ms_in_opt))
        print filename, 'opt', np.max(ms_in_opt), np.max(bheb_in_opt)

        rgb_ir = [l for l  in lines if l.startswith('rgb ir')]
        rgb_data = zip(*[t.strip().split()[2:] for t in rgb_ir])
        rgb_data = np.array(rgb_data, dtype=float)
        eagb_in_rgb = rgb_data[5]/rgb_data[7]
        rheb_in_rgb = rgb_data[4]/rgb_data[7]
        ir_eagb_contam = np.append(ir_eagb_contam, np.max(eagb_in_rgb))
        ir_rheb_contam = np.append(ir_rheb_contam, np.max(rheb_in_rgb))
        #print filename, 'ir', np.max(eagb_in_rgb), np.max(rheb_in_rgb)

        ir =  [l for l  in lines if l.startswith('rgb ir') or l.startswith('agb ir')]
        data = zip(*[t.strip().split()[2:] for t in ir])
        data = np.array(data, dtype=float)
        ms_in_ir = data[0]/data[7]
        bheb_in_ir = data[3]/data[7]
        #print filename, 'ir', np.max(eagb_in_rgb), np.max(rheb_in_rgb)

        ir_bheb_contam = np.append(ir_bheb_contam, np.max(bheb_in_ir))
        ir_ms_contam = np.append(ir_ms_contam, np.max(ms_in_ir))

    print 'opt eagb, rheb', np.max(opt_eagb_contam), np.max(opt_rheb_contam)
    print 'ir eagb, rheb', np.max(ir_eagb_contam), np.max(ir_rheb_contam)
    print 'opt bheb, ms', np.max(opt_bheb_contam), np.max(opt_ms_contam)
    print 'ir bheb, ms', np.max(ir_bheb_contam), np.max(ir_ms_contam)


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

    def poission_chi2(self, hist_it_up=False, table_file='default',
                      just_gauss=False):
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
        opt_model_hists, opt_models_binss, opt_model_norms = self.files.load_lf_file(self.fnames[0])
        ir_model_hists, ir_models_binss, ir_model_norms = self.files.load_lf_file(self.fnames[1])

        opt_chi2 = np.array([])
        ir_chi2 = np.array([])

        opt_chi2_agb = np.array([])
        ir_chi2_agb = np.array([])
        opt_pval = np.array([])
        ir_pval = np.array([])

        opt_pval_agb = np.array([])
        ir_pval_agb = np.array([])

        nhists = np.min([len(ir_model_hists), len(opt_model_hists)])

        for i in range(nhists):
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(opt_gal.hist[obins[1:]],
                                                           opt_model_hists[i][obins[1:]])
            if just_gauss is True:
                chi2, pval = self.chi2(opt_gal.hist[obins[1:]],
                                       opt_model_hists[i][obins[1:]])
                opt_pval = np.append(opt_pval, pval)

            opt_chi2 = np.append(opt_chi2, chi2)
            #print 'opt', chi2, np.mean(sig)
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(ir_gal.hist[ibins[1:]],
                                                           ir_model_hists[i][ibins[1:]])
            if just_gauss is True:
                chi2, pval = self.chi2(ir_gal.hist[ibins[1:]],
                                       ir_model_hists[i][ibins[1:]])
                ir_pval = np.append(ir_pval, pval)

            ir_chi2 = np.append(ir_chi2, chi2)
            #print 'ir', chi2, np.mean(sig)

            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(opt_gal.hist[agb_obins[1:]],
                                                           opt_model_hists[i][agb_obins[1:]])
            if just_gauss is True:
                chi2, pval = self.chi2(opt_gal.hist[agb_obins[1:]],
                                       opt_model_hists[i][agb_obins[1:]])
                opt_pval_agb = np.append(opt_pval_agb, pval)

            opt_chi2_agb = np.append(opt_chi2_agb, chi2)
            #print 'opta', chi2, np.mean(sig)

            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(ir_gal.hist[agb_ibins[1:]],
                                                           ir_model_hists[i][agb_ibins[1:]])
            if just_gauss is True:
                chi2, pval = self.chi2(ir_gal.hist[agb_ibins[1:]],
                                       ir_model_hists[i][agb_ibins[1:]])
                ir_pval_agb = np.append(ir_pval_agb, pval)
            ir_chi2_agb = np.append(ir_chi2_agb, chi2)
            #print 'ira ', chi2, np.mean(sig)
        if just_gauss is True:
            return opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb, opt_pval, ir_pval, opt_pval_agb, ir_pval_agb
        else:
            return opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb

    def chi2(self, ohist, mhist):
        from scipy import stats
        # t-test:
        #ti = (ohist - mhist) ** 2 / (ohist + mhist)
        #naners = np.isnan(ti)
        #ti[naners] = 0
        #print np.sum(ti)
        # maybe there is a better way to mask this... chiw was a typo...
        chiw = (mhist - ohist) ** 2 / ohist
        naners = np.isinf(chiw)
        chiw[naners] = 0
        naners = np.isnan(chiw)
        chiw[naners] = 0
        oinds, = np.nonzero(chiw > 0)
        chi2 = np.sum(chiw)/float(len(oinds)-1)#, len(oinds)
        pval = 1 - stats.chi2.cdf(chi2, len(oinds)-1)
        return chi2, pval


def result2dict(result_files, search=None):
    res_dict = {}

    if 'narratio' in result_files[0]:
        if search is None:
            search = 'ratio'
        chi2 = False
    elif 'chi2' in result_files[0]:
        if search is None:
            search = 'chi2'
        chi2 = True
    else:
        print 'either narratio file or chi2 file'
        return {}

    for result_file in result_files:
        data = rsp.fileIO.readfile(result_file, string_column=0)
        target = os.path.split(result_file)[1].split('_')[3]
        agb_mod = os.path.split(result_file)[1].split(target)[0][:-1]
        if search == 'all':
            fields = [d for d in data.dtype.names if not 'target' in d]
        else:
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
                     outfile_loc='default', extra_str='', just_gauss=False):
    if just_gauss is True:
        extra_str2 = extra_str + '_gauss'
    else:
        extra_str2 = extra_str
    chi2_files = []
    for target in targets:
        for cmd_input_file in cmd_input_files:
            st = StatisticalComparisons(cmd_input_file, target,
                                        outfile_loc=outfile_loc,
                                        extra_str=extra_str)
            chi2_file = os.path.join(outfile_loc,
                                     '%s_%s%s_chi2.dat' % (st.agb_mod, target,
                                                           extra_str2))
            result =  st.poission_chi2(table_file=table_file, just_gauss=just_gauss)
            if just_gauss is True:
                opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb, opt_pval, ir_pval, opt_pval_agb, ir_pval_agb = result
                cfmt = '%i %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n'
                with open(chi2_file, 'w') as c2:
                    c2.write('# sfr opt_chi2 ir_chi2 opt_agb_chi2 ir_agb_chi2 ')
                    c2.write('ir_chi2_agb opt_pval ir_pval opt_pval_agb ir_pval_agb \n')
                    for i in range(len(opt_chi2)):
                        c2.write(cfmt % (i, opt_chi2[i], ir_chi2[i],
                                         opt_chi2_agb[i], ir_chi2_agb[i],
                                         opt_pval[i], ir_pval[i],
                                         opt_pval_agb[i], ir_pval_agb[i]))
                chi2_files.append(chi2_file)
            else:
                opt_chi2, ir_chi2, opt_chi2_agb, ir_chi2_agb = result
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
    ir_table = np.empty(shape=(len(targets) + 1, len(agb_mods)*2 + 1),
                        dtype='|S20')
    ir_table[:, :] = ''
    opt_table = np.empty(shape=(len(targets) + 1, len(agb_mods)*2 + 1),
                         dtype='|S20')
    opt_table[:, :] = ''

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
        # column 0 is target, columns 1, 3, 5 have ratios
        column = (agb_mods.index(agb_mod) * 2) + 1
        #print row, column
        # target
        table[row, 0] = '%s &' % target.upper()

        # data
        dnarr = data_dict['%s_%s' % (target.lower(), band)]
        derr = data_dict['%s_%s_err' % (target.lower(), band)]
        #dstr = fmt % (dnarr, derr)
        #table[row, 1] = dstr

        # model
        mnarr = val
        # grab the error from the dict
        err_key = key.replace('mean', 'err_mean')
        mnerr = nar_dict[err_key]
        mstr = fmt % (mnarr, mnerr)
        table[row, column] = mstr

        # frac difference
        #pct_diff = (mnarr - dnarr) / dnarr
        #pct_diff_err = np.abs(pct_diff * (mnerr/mnarr + derr/dnarr))
        pct_diff = (mnarr / dnarr)
        pct_diff_err = np.abs(pct_diff * (mnerr/mnarr + derr/dnarr))

        f = fmt
        # if final column, put \\ not &
        if column + 1 == table.shape[1] - 1:
            f = fmt2
        pdstr = f % (pct_diff, pct_diff_err)
        table[row, column + 1] = pdstr

    # totals:
    nar_dict = result2dict(narratio_files, search='all')
    data_total = data_table(targets)
    data_total = np.array(data_total.translate(None, '\\$Total&').replace('pm', ' ').split(), dtype=float)

    for i, agb_mod in enumerate(agb_mods):
        for band, table in zip(['ir', 'opt'], [ir_table, opt_table]):
            f = fmt
            column = (i * 2) + 1
            nrgb = np.sum([v for k, v in nar_dict.items()
                           if agb_mod in k and band in k and 'rgb' in k])
            nagb = np.sum([v for k, v in nar_dict.items()
                           if agb_mod in k and band in k and 'agb' in k])
            ratio = nagb / nrgb
            err = galaxy_tests.count_uncert_ratio(nagb, nrgb)
            table[-1, column] = f % (ratio, err)
            column += 1
            if column == table.shape[1] - 1:
                f = fmt2
            # frac difference
            if band == 'ir':
                j = 6
            if band == 'opt':
                j = 2
            dratio = data_total[j]
            derr = data_total[j+1]
            pct_diff = (ratio / dratio)
            pct_diff_err = np.abs(pct_diff * (err/ratio + derr/dratio))
            #table[-1, column] = f % (pct_diff, pct_diff_err)
            #table[-1, 0] = 'Total & '
    # mean

    for i, agb_mod in enumerate(agb_mods):
        for table in [ir_table, opt_table]:
            f = fmt
            column = (i * 2) + 1
            val, err = \
                np.mean(np.array([l.translate(None, ' $&\\').split('pm')
                    for l in table[:, column][:-1]], dtype=float), axis=0)
            print agb_mod, err/val
            table[-1, column] = f % (val, err)
            column += 1
            if column == table.shape[1] - 1:
                f = fmt2
            val, err = \
                np.mean(np.array([l.translate(None, ' $&\\').split('pm')
                    for l in table[:, column][:-1]], dtype=float), axis=0)
            table[-1, column] = f % (val, err)
            table[-1, 0] = 'Mean & '

    # write the file
    ratio_str = '$\\frac{N_{\\rm TP-AGB}}{N_{\\rm RGB}}$'
    header = 'Target & '
    for i in range(len(agb_mods)):
        header += '%s %s & Frac. Difference & ' % (ratio_str,
                                                   agb_mods[i].split('_')[-1])
    header = header[:-2] + '\\\\ \n \\hline \n'

    outfile_dir = os.path.split(narratio_files[0])[0]
    for band, table in zip(['opt', 'ir'], [opt_table, ir_table]):
        outfile = os.path.join(outfile_dir, '%s_narratio_table.tex' % band)
        with open(outfile, 'w') as out:
            out.write(header)
            np.savetxt(out, table, fmt='%s')
    return ir_table, opt_table


def data_table(targets, table_file='default'):
    """latex data table """
    # target, av, dist, [opt: frac comp, trgb, nrgb, nagb ratio] [ir: ..]
    if table_file == 'default':
        table_file = snap_src + '/tables/ancients_0.1_0.2_galaxies.dat'
    ags = sfh_tests_multi_proc.AncientGalaxies()
    ags.read_trgb_table(table_file)
    data_dict = get_data(table_file=table_file)
    comp_data = tables.read_completeness_table()
    row = ''
    row2 = ''
    ts = list(set([t.lower() for t in targets]) & set(ags.data.target))
    assert len(ts) == len(targets), 'cant find all targets in ags.data'
    inds = list([np.where(t.lower() == ags.data.target)[0][0] for t in targets])
    opt_agb_tot = np.sum(ags.data.nopt_agb[inds])
    opt_rgb_tot = np.sum(ags.data.nopt_rgb[inds])
    ir_agb_tot = np.sum(ags.data.nir_agb[inds])
    ir_rgb_tot = np.sum(ags.data.nir_rgb[inds])

    totfmt = 'Total & %i & %i & $%.3f\\pm%.3f$ &  %i & %i & $%.3f\\pm%.3f$ \\\\'
    total = totfmt % (opt_agb_tot, opt_rgb_tot, opt_agb_tot/opt_rgb_tot,
                      galaxy_tests.count_uncert_ratio(opt_agb_tot, opt_rgb_tot),
                      ir_agb_tot, ir_rgb_tot, ir_agb_tot/ir_rgb_tot,
                      galaxy_tests.count_uncert_ratio(ir_agb_tot, ir_rgb_tot))

    for target in targets:
        if target.upper() == 'NGC2976-DEEP':
            extra_key='F606W,F814W'
        else:
            extra_key=None
        (Av, dmod) = [angst_data.get_item(target, i, extra_key=extra_key)
                      for i in ['Av', 'dmod']]
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
        row2 += '$%.3f\\pm%.3f$ & ' % (sub_dict['%s_opt' % target.lower()],
                                       sub_dict['%s_opt_err' % target.lower()])
        row2 += '%(nir_agb)i & %(nir_rgb)i & ' % nstars_row
        row2 += '$%.3f\\pm%.3f$ \\\\ \n' % (sub_dict['%s_ir' % target.lower()],
                                            sub_dict['%s_ir_err' % target.lower()])

    with open('data_table.tex', 'w') as out:
        out.write(row)
        out.write(row2)
        out.write(total)

        out.write('# max err pct. opt: %.3f ir: %.3f \n' % \
                  (np.max(opt_err_pct), np.max(ir_err_pct)))
    return total


def chi2plot(model_dict, outfile_loc=None):
    targets = list(np.unique([k.split('_')[0] for k in model_dict.keys()]))
    agb_mods = list(np.unique(['_'.join(k.split('_')[1:4])
                               for k in model_dict.keys()]))

    cols = ['darkgreen', 'navy', 'darkred']

    fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=False,
                            figsize=(10,10))
    offsets = np.linspace(0, 1, len(targets))
    for key, val in model_dict.items():
        if 'std' in key:
            continue
        target = key.split('_')[0]
        errval = model_dict[key.replace('mean', 'std')]
        ioff = targets.index(target)
        agb_mod = '_'.join(key.split('_')[1:4])
        col = cols[agb_mods.index(agb_mod)]
        sym = 'o'
        if not agb_mod.endswith('nov13'):
            mfc='white'
        else:
            mfc = col
        if not 'nov13' in agb_mod:
            sym = '*'
        ax_row = 0
        ax_col = 0
        if not '_agb' in key:
            ax_row = 1
        if 'ir' in key:
            ax_col = 1
        ax = axs[ax_row][ax_col]

        ax.errorbar(offsets[ioff], val, yerr=errval, marker=sym, color=col, ms=12,
                    mfc=mfc, ecolor='black', mew=1.5, elinewidth=2)
        ax.set_ylabel('$\chi^2$', fontsize=20)

        ax.xaxis.set_ticks(offsets)
        ax.set_xticklabels(['$%s$' % t.replace('-deep', '').replace('-', '\!-\!').upper() for t in targets])
        [t.set_rotation(30) for t in ax.get_xticklabels()]

        #ymaxs[ax_num] = np.max([val, ymaxs[ax_num]])
    axs[0][0].set_title(r'$\rm{Optical}$', fontsize=20)
    axs[0][1].set_title(r'$\rm{NIR}$', fontsize=20)
    [ax.set_ylim(0, 25) for ax in axs[:, 0]]
    [ax.set_ylim(0, 10) for ax in axs[:, 1]]
    fig.subplots_adjust(hspace=0.1)
    xlims = ax.get_xlim()
    off = np.diff(offsets)[0]
    ax.set_xlim(xlims[0]-off/2, xlims[1]+off/2)
    sym = ['o', 'o', '*']
    mfc = [cols[0], 'None', 'None']
    [axs[0, 0].plot(-99, 99, sym[j], mfc=mfc[j], ms=12, mew=1.5, color=cols[j],
     label='$%s$' % model_plots.translate_model_name(agb_mods[j].split('_')[-1]))
     for j in range(len(agb_mods))]
    axs[0, 0].legend(loc=0, numpoints=1)
    [ax.annotate(r'$\rm{TP\!-\!AGB\ Only}$', (0.02, 0.02), fontsize=16,
                 xycoords='axes fraction') for ax in axs[0, :]]
    if outfile_loc is None:
        outfile_loc = os.getcwd()
    outfile = os.path.join(outfile_loc, 'chi2_plot.png')
    plt.tick_params(labelsize=16)
    plt.savefig(outfile, dpi=150)

    return axs


def run_match_stats(targets='ancients'):
    targets = galaxy_tests.load_targets(targets)
    hmc_file_loc = os.path.join(snap_src, 'data', 'sfh_parsec')
    cmd_file_loc = os.path.join(hmc_file_loc, 'cmd_files')
    for target in targets:
        target = target.lower()
        try:
            target = target.replace('-deep', '')
            hmc_file, = rsp.fileIO.get_files(hmc_file_loc, '%s*sfh' % target)
            cmd_file, = rsp.fileIO.get_files(cmd_file_loc, '%s*cmd' % target)
        except:
            print target, 'cmd file not found.'
            continue
        rsp.match_utils.match_stats(hmc_file, cmd_file)
    return