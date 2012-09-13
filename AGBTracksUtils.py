import pdb
import trilegal_diagnostics
import GoogleSitesTable
import galaxy_tests
import os
import sys
import numpy as np
import fileIO
import graphics
import matplotlib.pyplot as plt


def get_TP_inds(ntp):
    un = np.unique(ntp)
    # this is the first step in each TP.
    iTPs = [list(ntp).index(u) for u in un]
    # The indices of each TP.
    TPs = [np.arange(iTPs[i], iTPs[i + 1]) for i in range(len(iTPs) - 1)]
    # don't forget the last one.
    TPs.append(np.arange(iTPs[i + 1], len(ntp)))
    return TPs


def get_unique_inds(ntp):
    un = np.unique(ntp)
    # this is the first step in each TP.
    iTPs = [list(ntp).index(u) for u in un]
    # The indices of each TP.
    TPs = [np.arange(iTPs[i], iTPs[i + 1]) for i in range(len(iTPs) - 1)]
    # don't forget the last one.
    TPs.append(np.arange(iTPs[i + 1], len(ntp)))
    return TPs, iTPs


def add_points_to_q_track(track, qs):
    '''
    when to add an extra point for low masses
    if logt[qs+1] is hotter than logt[qs]
    and there is a point inbetween logt[qs] and logt[qs+1] that is cooler
    than logt[qs]
    add the coolest point.
    '''
    addpt = []
    logt = track.get_col('T_star')
    step = track.get_col('step')
    status = track.get_col('status')
    Tqs = logt[qs]
    # need to use some unique array, not log t, since log t could repeat,
    # index would find the first one, not necessarily the correct one.
    Sqs = step[qs] - 1.  # makes it the same as qs.
    # takes the difference in logt(qs) to see if we get hotter or colder.
    # finds where the logt goes from getting colder to hotter...
    ht = np.where(np.diff(np.sign(np.diff(Tqs))))[0] + 1
    ht = np.append(ht, ht + 1)  # between the first and second
    Sqs_ht = Sqs[ht]
    # the indices between each hot point.
    t_mids = [map(int, step[int(Sqs_ht[i]): int(Sqs_ht[i + 1])])
              for i in range(len(Sqs_ht) - 1)]
    Sqs_ht = Sqs_ht[: -1]
    for i in range(len(Sqs_ht)):
        hot_inds = np.nonzero(logt[int(Sqs_ht[i])] > logt[t_mids[i]])[0]
        if len(hot_inds) > 0:
            # index of the min T of the hot index from above.
            addpt.append(list(logt).index(np.min(logt[[t_mids[i][hi]
                                                       for hi in hot_inds]])))

    if len(addpt) > 0:
        addpt = np.unique([a for a in addpt if status[a] == 7.])
    Qs = np.sort(np.concatenate((addpt, qs)))
    return Qs, addpt


def find_dldt(track, TPs, addpt):
    '''
    Finds dL/dt of track object
    '''
    status = track.get_col('status')
    #dl/dT seems somewhat linear for 0.2 < phi < 0.4 ...
    phi = track.get_col('PHI_TP')
    # The first line in the agb track is 1.This isn't the quiessent...
    phi[0] = -1.
    lin_rise, = np.nonzero((status == 7) & (phi < 0.4) & (phi > 0.2))
    rising = [list(set(TP) & set(lin_rise)) for TP in TPs]
    logl = fileIO.AGBTracks.get_col(track, 'L_star')
    logt = fileIO.AGBTracks.get_col(track, 'T_star')
    order = 1
    fits = [np.polyfit(logt[r], logl[r], order) for r in rising if len(r) > 0]
    slopes = [fits[i][0] for i in range(len(fits))]
    if len(addpt) > 0:
        addrise = [TPs.index(TP) for TP in TPs
                   if len(set(addpt) & set(TP)) > 0]
        addslope = slopes[addrise[0]]
        Slopes = []
        for s in slopes:
            if s == addslope:
                Slopes.append(addslope)
            Slopes.append(s)
    else:
        Slopes = slopes[:]
    return rising, Slopes, fits


def do_everything(**kwargs):
    '''
    HOW TO USE THIS:
    This script formats Paola's tracks and creates the files needed to use them
    with TRILEGAL.

    You can either pass a dictionary with the necessary locations of the files
    to be made and the files made or simply make an input file:

    kwargs = load_input(filename)
    do_everything(**kwargs)

    the input file should be formatted in this way (comments are optional):
    # Paola's tracks are here:
    agbtrack_dir    /Users/phil/research/TP-AGBcalib/AGBTracks
    # This agb mix name
    agb_mix         agb_caf09
    # This set name
    set_name        S0
    # Where cmd_input files will go for trilegal
    trilegal_dir    output_for_trilegal/
    # Paola's formatted tracks are here:
    isotrack_dir    isotrack/parsec/AGB_TRACKS/
    # File to link from cmd_input to Paola's formatted tracks:
    tracce_dir      isotrack/parsec/
    # make initial and final mass relation? (Need more than one metallicity)
    make_imfr       False
    # Diagnostic plots will mimic directory structure of abgtrack_dir
    diagnostic_dir  DIAGPLOTS

    WHAT'S GOING ON:
    Here is what it takes to go from Paola's agb tracks to trilegal:

    1. Paola's tracks
    2. Paola's tracks formatted for trilegal
    3. trilegal files that link to Paola's formatted tracks

    1. Paola's tracks should have naming scheme
        [mix]/[set]/[Trash]_[metallicity]_[Y]/[track_identifier]
    ex agb_caf09_z0.008/S1/agb_*Z*.dat
    ex CAF09/S_SCS/*Z*/agb_*Z*.dat
    These can go anywhere, and do not need to be included in trilegal
    directories.

    2. Paola's tracks are formatted for trilegal with no header file.
    These go in trilegal_1.3/isotrack/AGB_TRACKS/Z[metallicity]_[mix]_[set].dat
    They only include the quiescent phases and also a column with dL/dT

    3. trilegal needs two files to link to Paola's formatted tracks
       a. track file goes here: trilegal_1.3/isotrack/tracce_[mix]_[set].dat
       b. cmd_input_file that links to the track file
          goes here trilegal_1.3/cmd_input_[mix]_[set].dat
    '''
    home = os.getcwd()
    over_write = kwargs.get('overwrite', False)
    agbtrack_dir = kwargs.get('agbtrack_dir')
    agb_mix = kwargs.get('agb_mix')
    set_name = kwargs.get('set_name')
    trilegal_dir = kwargs.get('trilegal_dir')
    isotrack_dir = kwargs.get('isotrack_dir')
    tracce_dir = kwargs.get('tracce_dir')
    diagnostic_dir0 = kwargs.get('diagnostic_dir')
    make_imfr = kwargs.get('make_imfr')
    make_copy = kwargs.get('make_copy')
    mass_loss = kwargs.get('mass_loss')
    # do dirs exist?
    fileIO.ensure_dir(isotrack_dir)
    fileIO.ensure_dir(os.path.join(diagnostic_dir0, agb_mix, set_name + '/'))

    # name convention: [mix]_[set].dat
    name_conv = '_'.join((agb_mix, set_name)) + '.dat'

    if make_imfr is True:
        ifmr = os.path.join(diagnostic_dir0, agb_mix, set_name,
                            '_'.join(('ifmr', name_conv)))
        mfile = open(ifmr, 'w')
        mfile.write('# M_i M_f Z\n')

    # cmd_input_file that links to the track file
    cmd_input = '_'.join(('cmd', 'input', name_conv))
    cmd_input_file = os.path.join(home, trilegal_dir, cmd_input)

    # track file to link from cmd_input to paola's formatted tracks
    tracce_fh = '_'.join(('tracce', name_conv))
    tracce_file = os.path.join(home, tracce_dir, tracce_fh)
    tracce_file_rel = os.path.join('isotrack', tracce_dir, tracce_fh)
    # to find AGB tracks (this could be a kwarg...
    # I have Z because agb_imfr.dat would crash.
    track_identifier = 'agb_*Z*.dat'

    isofiles, Zs, Ys = [], [], []

    # moving to the the directory with metallicities.
    working_dir = os.path.join(agbtrack_dir, agb_mix, set_name)
    os.chdir(working_dir)
    metal_dirs = [m for m in os.listdir(working_dir) if os.path.isdir(m)]
    metals = np.argsort([float(m.split('_')[1].replace('Z', ''))
                         for m in metal_dirs])
    print 'found %i metallicities' % len(metal_dirs)

    lifetimesfile = os.path.join(diagnostic_dir0, agb_mix, set_name,
                                 '_'.join(('tau_cm', name_conv)))

    cm = open(lifetimesfile, 'w')
    cm.write('# z mass tauc taum\n')
    for metal_dir in np.array(metal_dirs)[metals]:
        if diagnostic_dir0 is not None:
            diagnostic_dir = os.path.join(diagnostic_dir0, agb_mix, set_name,
                                          metal_dir) + '/'
            fileIO.ensure_dir(diagnostic_dir)
            # update kwarg
            kwargs['diagnostic_dir'] = diagnostic_dir

        agb_tracks = fileIO.get_files(os.path.join(working_dir, metal_dir),
                               track_identifier)
        agb_tracks.sort()
        ax = plt.axes()
        kwargs['ax'] = ax
        # Paola's formatted tracks
        metallicity = float(metal_dir.split('_')[1].replace('Z', ''))
        Y = float(metal_dir.split('_')[-1].replace('Y', ''))
        iso_name_conv = '_'.join(('Z%.4f' % metallicity, name_conv))
        isofile = os.path.join(home, isotrack_dir, iso_name_conv)
        if over_write is False and os.path.isfile(isofile):
            print 'not over writing %s' % isofile
            continue
        isofile_rel_name = os.path.join('isotrack', isotrack_dir,
                                        iso_name_conv)
        out = open(isofile, 'w')
        out.write(('# age(yr) logL logTe m_act mcore c/o period ip Mdot(Msun/yr) logTe X Y CNO dlogTe/dlogL \n'))
        print 'found %i tracks' % len(agb_tracks)
        for agb_track in agb_tracks:
            flag = 0
            addpt = []
            # load track
            track = fileIO.get_numeric_data(agb_track)
            tdata = track.data_array
            ntp = tdata['NTP']
            age = tdata['ageyr']
            phi = tdata['PHI_TP']
            logt = tdata['T_star']
            logl = tdata['L_star']
            dt = tdata['dt']
            co = tdata['CO']
            mdot = tdata['dMdt']
            mass = track.mass
            met = metallicity
            mstar, = np.nonzero((co <= 1) & (logl >= 3.3) & (mdot <= -5))
            cstar, = np.nonzero((co >= 1) & (mdot <= -5))
            try:
                tauc = np.sum(dt[cstar]) / 1e6
            except IndexError:
                tauc = 0.
            try:
                taum = np.sum(dt[mstar]) / 1e6
            except IndexError:
                taum = 0.
            cm.write(' %.4f %.3f %.4f %.4f\n' % (met, mass, tauc, taum))
            # get the indices of the thermal pulses.
            if len(ntp) == sum(ntp):
                print 'no tracks!', agb_track
                flag = 1
            if flag == 1:
                continue
            TPs = get_TP_inds(ntp)
            # The first line in the agb track is 1. This isn't the quiescent...
            phi[0] = -1.
            # The quiescent phase is the the max phase in each TP,
            # i.e., closest to 1.

            qs = np.cumsum([np.argmax(phi[TP]) for TP in TPs])
            Qs = qs
            # when to add an extra point to quiescent track
            if len(Qs) <= 9 and mass < 3.:
                Qs, addpt = add_points_to_q_track(track, qs)
            Qs = list(Qs)
            # find dl/dt of track
            rising, slopes, fits = find_dldt(track, TPs, addpt)

            # need to add first line.
            logt0 = logt[0]
            logl0 = logl[0]
            logt1 = logt[2]
            logl1 = logl[2]
            slope = (logl1 - logl0) / (logt1 - logt0)
            Qs.insert(0, 0)
            slopes.insert(0, slope)

            # make diagnostic plots and imfr
            if diagnostic_dir is not None:
                graphics.diag_plots(track, logl, logt, age, slopes, Qs, addpt,
                                    rising, fits, **kwargs)
            # make iso file for trilegal
            fileIO.make_iso_file(track, Qs, slopes, out)
            if kwargs['make_imfr'] is True:
                M_s = fileIO.AGBTracks.get_col(track, 'M_star')
                mfile.write(' %f %f %f\n' % (M_s[0], M_s[-1], float(met)))
        out.close()
        print 'wrote', isofile
        fileIO.make_local_copy(isofile, dest=make_copy)
        # keep information for tracce file
        isofiles.append(isofile_rel_name)
        Ys.append(Y)
        Zs.append(metallicity)

    # make file to link cmd_input to formatted agb tracks
    metfile = fileIO.make_met_file(tracce_file, Zs, Ys, isofiles)
    print 'wrote', tracce_file
    fileIO.make_local_copy(tracce_file, dest=make_copy)
    # make cmd_input file that
    cmd_input = fileIO.write_cmd_input_file(**{'cmd_input_file': cmd_input_file,
                                            'file_tpagb': tracce_file_rel,
                                            'mass_loss': mass_loss})
    print 'wrote', cmd_input_file
    cm.close()
    print 'wrote', lifetimesfile
    graphics.plot_cluster_test(lifetimesfile, **kwargs)
    if make_imfr is True:
        mfile.close()
        print 'wrote', ifmr
        graphics.plot_ifmr(ifmr)

    os.chdir(home)
    return cmd_input_file


if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
    except:
        print do_everything.__doc__
    pdb.set_trace()
    kwargs = fileIO.load_input(input_file)
    # Paola's tracks -> trilegal + tests based only on tracks
    overwrite = kwargs.get('overwrite')
    kwargs['overwrite'] = True
    cmd_input_file = do_everything(**kwargs)

    agb_mix = kwargs['agb_mix']
    set_name = kwargs['set_name']
    track_set = '_'.join((agb_mix, set_name))

    diagnostic_dir = kwargs.get('diagnostic_dir')
    # Marco's scripts to run trilegal at age and z

    if kwargs.get('trilegal_diagnostics'):
        tri_dir = os.path.join(diagnostic_dir, 'trilegal_files')
        sfh_dir = os.path.join(tri_dir, 'sfh')
        plt_dir = os.path.join(diagnostic_dir, agb_mix, set_name, track_set)
        trilegal_diagnostics.main(track_set, sfh_dir, tri_dir, plt_dir,
                                  over_write=kwargs['overwrite'])
        if kwargs.get('google_table'):
            image_location = kwargs.get('image_location')
            if not image_location:
                image_location = plt_dir
            GoogleSitesTable.trilegal_diag_table(image_location)

    # Scripts to make LF compared to data
    if kwargs.get('IDs'):
        gt_kw = {'outdir': kwargs['galaxy_outdir']}
        if kwargs.get('google_table'):
           gt_kw['make_plots'] = True
           gt_kw['publish_plots'] = True
        
        galaxy_tests.main(kwargs['IDs'], ['%s.dat' % track_set], **gt_kw)
        pass
