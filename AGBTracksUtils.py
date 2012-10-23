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
import multiprocessing



def metallicity_from_dir(met):
    z = float(met.split('_')[1].replace('Z', ''))
    y = float(met.split('_')[-1].replace('Y', ''))
    return z, y


def AGB_file_setup(infile):
    
    # Check for Paola's formatted tracks
    fileIO.ensure_dir(infile.isotrack_dir)

    # are we making diagnostic plots, check directory.
    if infile.diagnostic_dir0:
        fileIO.ensure_dir(os.path.join(infile.diagnostic_dir0,
                                       infile.agb_mix,
                                       infile.set_name + '/'))
    else:
        print 'not making diagnostic plots'

    # set name convention: [mix]_[set].dat
    infile.name_conv = '_'.join((infile.agb_mix, infile.set_name)) + '.dat'

    # set track search string
    infile.track_identifier = 'agb_*Z*.dat'

    # cmd_input_file that links to the track file
    cmd_input = '_'.join(('cmd', 'input', infile.name_conv))
    infile.cmd_input_file = os.path.join(infile.home, infile.trilegal_dir, cmd_input)

    # track file to link from cmd_input to paola's formatted tracks
    tracce_fh = '_'.join(('tracce', infile.name_conv))
    infile.tracce_file = os.path.join(infile.home, infile.tracce_dir, tracce_fh)
    infile.tracce_file_rel = os.path.join('isotrack', infile.tracce_dir, tracce_fh)

    # moving to the the directory with metallicities.
    infile.working_dir = os.path.join(infile.agbtrack_dir, infile.agb_mix, infile.set_name)
    os.chdir(infile.working_dir)
    metal_dirs = [m for m in os.listdir(infile.working_dir) if os.path.isdir(m)]
    metals = np.argsort([metallicity_from_dir(m)[0] for m in metal_dirs])
    infile.metal_dirs = np.array(metal_dirs)[metals]
    print 'found %i metallicities' % len(metal_dirs)


def multipro_plots(agb_tracks, infile):
    for agb_track in agb_tracks:
        # load track
        track = fileIO.get_numeric_data(agb_track)
        if track.bad_track is True:
            continue
        # make diagnostic plots
        graphics.diag_plots(track, infile)    
    return
    
def make_plots(infile):
    if infile.diagnostic_dir0 is None:
        return 

    pool = multiprocessing.Pool()
    res = []

    for metal_dir in infile.metal_dirs:
        diagnostic_dir = os.path.join(infile.diagnostic_dir0,
                                      infile.agb_mix,
                                      infile.set_name,
                                      metal_dir) + '/'
        fileIO.ensure_dir(diagnostic_dir)
        # update infile class to place plots in this directory
        infile.diagnostic_dir = diagnostic_dir
        agb_tracks = fileIO.get_files(os.path.join(infile.working_dir, metal_dir),
                                      infile.track_identifier)
        agb_tracks.sort()

        res.append(pool.apply_async(multipro_plots, (agb_tracks, infile)))

    for r in res:
        r.get()

    return

def do_everything(infile):
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
    infile.home = os.getcwd()

    # set up file names and directories, cd to paola's tracks.
    AGB_file_setup(infile)

    isofiles, Zs, Ys = [], [], []
    for metal_dir in infile.metal_dirs:
        metallicity, Y = metallicity_from_dir(metal_dir)
        if infile.diagnostic_dir0 is not None:
            diagnostic_dir = os.path.join(infile.diagnostic_dir0,
                                          infile.agb_mix,
                                          infile.set_name,
                                          metal_dir) + '/'
            fileIO.ensure_dir(diagnostic_dir)
            # update infile class to place plots in this directory
            infile.diagnostic_dir = diagnostic_dir

        agb_tracks = fileIO.get_files(os.path.join(infile.working_dir, metal_dir),
                                      infile.track_identifier)
        agb_tracks.sort()

        iso_name_conv = '_'.join(('Z%.4f' % metallicity, infile.name_conv))
        isofile = os.path.join(infile.home, infile.isotrack_dir, iso_name_conv)

        if infile.over_write is False and os.path.isfile(isofile):
            print 'not over writing %s' % isofile
            continue

        isofile_rel_name = os.path.join('isotrack', infile.isotrack_dir,
                                        iso_name_conv)
        imfr_data = np.array([])
        lifetime_data = np.array([])
        with open(isofile, 'w') as out:
            out.write('# age(yr) logL logTe m_act mcore c/o period ip')
            out.write(' Mdot(Msun/yr) logTe X Y CNO dlogTe/dlogL \n')
            print 'found %i tracks' % len(agb_tracks)
            for agb_track in agb_tracks:
                # load track
                track = fileIO.get_numeric_data(agb_track)
                if track.bad_track is True:
                    continue

                # make iso file for trilegal
                fileIO.make_iso_file(track, out)
                # save information for lifetime file.
                lifetime_datum = np.array([metallicity, track.mass, track.tauc,
                                           track.taum])
                lifetime_data = np.append(lifetime_data, lifetime_datum)

                # make diagnostic plots
                if infile.diagnostic_dir0 is not None:
                    graphics.diag_plots(track, infile)
    
                # save information for imfr
                if infile.make_imfr is True:
                    M_s = track.get_col('M_star')
                    imfr_datum = np.array([M_s[0], M_s[-1],
                                           float(metallicity)])
                    imfr_data = np.append(imfr_data, imfr_datum)

        print 'wrote', isofile
        fileIO.make_local_copy(isofile, dest=infile.make_copy)
        # keep information for tracce file
        isofiles.append(isofile_rel_name)
        Ys.append(Y)
        Zs.append(metallicity)

    # make file to link cmd_input to formatted agb tracks
    metfile = fileIO.make_met_file(infile.tracce_file, Zs, Ys, isofiles)

    fileIO.make_local_copy(infile.tracce_file, dest=infile.make_copy)

    # make cmd_input file
    cmd_input = fileIO.write_cmd_input_file(**{'cmd_input_file': infile.cmd_input_file,
                                               'file_tpagb': infile.tracce_file_rel,
                                               'mass_loss': infile.mass_loss})

    make_plots(infile)

    if infile.make_imfr is True and infile.diagnostic_dir0 is not None:
        ifmr_file = os.path.join(infile.diagnostic_dir0, infile.agb_mix, infile.set_name,
                                 '_'.join(('ifmr', name_conv)))
        fileIO.savetxt(ifmr_file, imfr_data, header='# M_i M_f Z\n')
        graphics.plot_ifmr(ifmr_file)

    if infile.diagnostic_dir0 is not None:
        lifetime_file = os.path.join(infile.diagnostic_dir0, infile.agb_mix, infile.set_name,
                                     '_'.join(('tau_cm', name_conv)))
        fileIO.savetxt(lifetime_file, lifetime_data, header='# z mass tauc taum\n')
        graphics.plot_cluster_test(lifetimesfile, **kwargs)

    os.chdir(infile.home)
    return cmd_input_file


if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
    except:
        print do_everything.__doc__
    pdb.set_trace()
    infile = fileIO.input_file(input_file)
    # Paola's tracks -> trilegal + tests based only on tracks
    cmd_input_file = do_everything(infile)

    agb_mix = infile.agb_mix
    set_name = infile.set_name
    track_set = '_'.join((agb_mix, set_name))

    diagnostic_dir = infile.diagnostic_dir0
    # Marco's scripts to run trilegal at age and z

    if infile.trilegal_diagnostics:
        tri_dir = os.path.join(diagnostic_dir, 'trilegal_files')
        sfh_dir = os.path.join(tri_dir, 'sfh')
        plt_dir = os.path.join(diagnostic_dir, agb_mix, set_name, track_set)
        trilegal_diagnostics.main(track_set, sfh_dir, tri_dir, plt_dir,
                                  over_write=infile.overwrite)
        if infile.google_table:
            image_location = infile.image_location
            if not image_location:
                image_location = plt_dir
            GoogleSitesTable.trilegal_diag_table(image_location)

    # Scripts to make LF compared to data
    if infile.IDs:
        gt_kw = {'outdir': infile.galaxy_outdir}
        if infile.google_table:
           gt_kw['make_plots'] = True
           gt_kw['publish_plots'] = True
        
        galaxy_tests.main(infile.IDs, ['%s.dat' % track_set], **gt_kw)
        pass
