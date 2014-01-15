import os
import numpy as np
import sfh_tests
import ResolvedStellarPops as rsp
from TPAGBparams import snap_src
import galaxy_tests

def ms_color_cut():
    comp90 = sfh_tests.read_completeness_table(absmag=True)
    tri_dir = os.environ['TRILEGAL_ROOT']
    # make these simulations by editing a constant sf trilegal file
    # using tab_sfr/ as a source of templates.
    sgals = [rsp.Galaxies.simgalaxy(tri_dir + 'const_sfr_out_z06.dat',
                                    filter1='F606W', filter2='F814W'),
             rsp.Galaxies.simgalaxy(tri_dir + 'const_sfr_out_z0006.dat',
                                    filter1='F606W', filter2='F814W')]
    dline = {}
    for band in ['opt', 'ir']:
        if band == 'opt':
            filter2 = 'F814W'
        else:
            filter2 = 'F160W'
        for sgal in sgals:
            sgal.all_stages()
            if not sgal.name in dline.keys():
                dline[sgal.name] = {}
            # only consider MS stars that are brighter than the faintest
            # 90% complteness mag in the sample
            iblue, = np.nonzero(sgal.data.get_col(filter2) <
                                np.max(comp90['%s_filter2' % band]))
            ims = list(set(iblue) & set(sgal.ims))
            ibheb = list(set(iblue) & set(sgal.ibheb))
            ims = list(np.concatenate([ims, ibheb]))
            if len(ims) == 0:
                # could be an issue that no MS stars are around...
                print 'warning', filter2, 'no MS found.'
                continue
            if band == 'opt':
                dline[sgal.name]['F606W'] = np.max(sgal.data.get_col('F606W')[ims])
                dline[sgal.name]['F814W'] = np.max(sgal.data.get_col('F814W')[ims])
                dline[sgal.name]['F475W'] = np.max(sgal.data.get_col('F475W')[ims])
            else:
                dline[sgal.name]['F110W'] = np.max(sgal.data.get_col('F110W')[ims])
                dline[sgal.name]['F160W'] = np.max(sgal.data.get_col('F160W')[ims])
    return dline

def color_cut_per_galaxy(table='default'):
    ms_dict = ms_color_cut()
    comp90 = sfh_tests.read_completeness_table(table=table, absmag=True)
    comp90_uncert = sfh_tests.read_completeness_table(table=table, uncertainties=True)
    #if table == 'default':
    #    table = snap_src + '/tables/completeness_0.90.dat'
    #table = table.replace('.dat', '_colorcuts.dat')
    photsys = 'wfc3snap'
    #fmt = '%(target)s %(opt_colorcut).3f %(ir_colorcut).3f \n'
    
    for sgalname, dline in ms_dict.items():
        print sgalname
        for i, target in enumerate(comp90['target']):
            for band in ['opt', 'ir']:
                if band == 'opt':
                    filter1 = sfh_tests.get_filter1(target.lower())
                    filter2 = 'F814W'
                else:
                    filter1 = 'F110W'
                    filter2 = 'F160W'
                m2m = {'target': target, 'filter2': filter2, 'filter1': filter1}
                if not filter1 in dline.keys():
                    continue
                Mag1 = dline[filter1]
                Mag2 = dline[filter2]
                mag1 = rsp.astronomy_utils.Mag2mag(Mag1, filter1, photsys, **m2m)
                mag2 = rsp.astronomy_utils.Mag2mag(Mag2, filter2, photsys, **m2m)
                color = mag1 - mag2
                color_uncert = comp90_uncert[i]['%s_color' % band]
                print target, filter1, filter2, '%.2f' % (color+color_uncert)
    
def find_contamination_by_phases(output_files=None):
    if output_files is None:
        output_files = [ snap_src + '/models/varysfh/ddo71/caf09_s_nov13/mc/output_ddo71_caf09_s_nov13.dat',
                        #snap_src + '/models/varysfh/ddo78/caf09_s_nov13/mc/output_ddo78_caf09_s_nov13.dat',
                        snap_src + '/models/varysfh/hs117/caf09_s_nov13/mc/output_hs117_caf09_s_nov13.dat',
                        #snap_src + '/models/varysfh/kdg73/caf09_s_nov13/mc/output_kdg73_caf09_s_nov13.dat',
                        snap_src + '/models/varysfh/kkh37/caf09_s_nov13/mc/output_kkh37_caf09_s_nov13.dat']#,
                        #snap_src + '/models/varysfh/ngc2976-deep/caf09_s_nov13/mc/output_ngc2976-deep_caf09_s_nov13.dat',
                        #snap_src + '/models/varysfh/ngc404/caf09_s_nov13/mc/output_ngc404_caf09_s_nov13.dat']

    for output_file in output_files:
        target = output_file.split('output_')[1].split('_')[0]
        print target
        filter1 = sfh_tests.get_filter1(target)

        ds = sfh_tests.Diagnostics(VarySFH_kw={'target': target})
        ds.mc = False

        sgal = rsp.Galaxies.simgalaxy(output_file, filter1=filter1,
                                      filter2='F814W')
        sgal.target = target

        sopt_rgb, sopt_agb, sir_rgb, sir_agb = \
            ds.do_normalization(filter1=filter1, trilegal_output=output_file,
                                hist_it_up=False, dry_run=True)

        ds.contamination_by_phases(sopt_rgb, sopt_agb, sir_rgb, sir_agb)

    return


def completeness_table_absmag(table='default'):
    '''
    convert the completeness table mags to abs mag.
    outfile = [table]_absmag.dat
    '''
    comp90 = sfh_tests.read_completeness_table(table)
    if table == 'default':
        table = snap_src + '/tables/completeness_0.90.dat'
    table = table.replace('.dat', '_absmag.dat')
    photsys = 'wfc3snap'
    fmt = '%(target)s %(opt_filter1).3f %(opt_filter2).3f %(ir_filter1).3f %(ir_filter2).3f \n'
    
    with open(table, 'w') as out:
        for i, target in enumerate(comp90['target']):
            dline = {'target': target}
            for band in ['opt', 'ir']:
                if band == 'opt':
                    filter1 = sfh_tests.get_filter1(target.lower())
                    filter2 = 'F814W'
                else:
                    filter1 = 'F110W'
                    filter2 = 'F160W'
                m2m = {'target': target, 'filter2': filter2, 'filter1': filter1}
                
                compf1 = comp90[i]['%s_filter1' % band]
                dline['%s_filter1' % band] = rsp.astronomy_utils.mag2Mag(compf1, filter1, photsys, **m2m)
                
                compf2 = comp90[i]['%s_filter2' % band]
                dline['%s_filter2' % band] = rsp.astronomy_utils.mag2Mag(compf2, filter2, photsys, **m2m)
            out.write(fmt % dline)
    print 'wrote %s' % table

def uncertainties_at_completeness(table='default', binwidth=0.1):
    '''
    write a table with the median uncertainties around the completeness value
    from the completeness table.
    '''
    comp90 = sfh_tests.read_completeness_table(table)
    if table == 'default':
        table = snap_src + '/tables/completeness_0.90.dat'
    table = table.replace('.dat', '_uncertainties.dat')
    opt_fits_src = snap_src + '/data/angst_no_trim'
    fmt = '%(target)s %(opt_filter1).3f %(opt_filter2).3f %(opt_color).3f %(ir_filter1).3f %(ir_filter2).3f %(ir_color).3f \n'
    title = '# ' + fmt.replace('%','').replace(')', '').replace('.3f','').replace('s','').replace('(','')
    with open(table, 'w') as out:
        out.write('# median uncertainty within +/-%.2f of completeness mag\n' % (binwidth/2))
        out.write(title)
        for i, target in enumerate(comp90['target']):
            ir_gal = galaxy_tests.load_galaxy(target, band='ir')
            opt_gal = galaxy_tests.load_galaxy(target, band='opt',
                                               fits_src=opt_fits_src)
            dline = {'target': target}
            for band, gal in zip(['opt', 'ir'], [opt_gal, ir_gal]):
                key = '%s_filter1' % band
                uncerts1, = np.nonzero((gal.mag1 < comp90[i][key] + binwidth/2) &
                                       (gal.mag1 > comp90[i][key] - binwidth/2))
                med_unct1 = np.median(gal.data.MAG1_ERR[uncerts1])
                dline[key] = med_unct1
                key = '%s_filter2' % band
                uncerts2, = np.nonzero((gal.mag2 < comp90[i][key] + binwidth/2) &
                                       (gal.mag2 > comp90[i][key] - binwidth/2))            
                
                med_unct2 = np.median(gal.data.MAG2_ERR[uncerts2])
                dline[key] = med_unct2
                med_color_unct = np.sqrt(med_unct1 ** 2 + med_unct2 ** 2)
                dline['%s_color' % band] = med_color_unct
                
            out.write(fmt % dline)
            
            
            
            
            
            
            
            
            
            
            
            