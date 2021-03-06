import pdb
import trilegal_diagnostics
import os
import sys
import numpy as np
import fileIO
import graphics
import ResolvedStellarPops as rsp

def metallicity_from_dir(met):
    ''' take Z and Y values from string'''
    if met.endswith('/'):
        met = met[:-1]

    if len(os.path.split(met)) > 0:
        met = os.path.split(met)[1]

    z = float(met.split('_')[1].replace('Z', ''))
    y = float(met.split('_')[-1].replace('Y', ''))

    return z, y


def AGB_file_setup(infile):
    '''set up files and directories for TPAGB parsing.'''
    infile.home = os.getcwd()

    # Check for Paola's formatted tracks
    rsp.fileIO.ensure_dir(infile.isotrack_dir)

    # are we making diagnostic plots, check directory.
    if infile.diagnostic_dir0:
        rsp.fileIO.ensure_dir(os.path.join(infile.diagnostic_dir0,
                                           infile.agb_mix,
                                           infile.set_name + '/'))
    else:
        print 'not making diagnostic plots'

    # set name convention: [mix]_[set].dat
    infile.name_conv = '%s.dat' % '_'.join((infile.agb_mix, infile.set_name))

    # set track search string
    infile.track_identifier = 'agb_*Z*.dat'

    # cmd_input_file that links to the track file
    cmd_input = '_'.join(('cmd', 'input', infile.name_conv))
    infile.cmd_input_file = os.path.join(infile.home, infile.trilegal_dir,
                                         cmd_input)

    # track file to link from cmd_input to paola's formatted tracks
    tracce_fh = '_'.join(('tracce', infile.name_conv))

    infile.tracce_file = os.path.join(infile.home, infile.tracce_dir,
                                      tracce_fh)

    infile.tracce_file_rel = os.path.join('isotrack', infile.tracce_dir,
                                          tracce_fh)

    # moving to the the directory with metallicities.
    infile.working_dir = os.path.join(infile.agbtrack_dir, infile.agb_mix,
                                      infile.set_name)
    os.chdir(infile.working_dir)
    metal_dirs = [m for m in os.listdir(infile.working_dir)
                  if os.path.isdir(m) and 'Z' in m]

    if infile.metals_subset is not None:
        print 'doing a subset of metallicities'
        metal_dirs = [m for m in metal_dirs
                      if metallicity_from_dir(m)[0] in infile.metals_subset]
    metals = np.argsort([metallicity_from_dir(m)[0] for m in metal_dirs])
    infile.metal_dirs = np.array(metal_dirs)[metals]
    print 'found %i metallicities' % len(metal_dirs)


def do_everything(infile):
    '''
    This script formats Paola's tracks and creates the files needed to use them
    with TRILEGAL as well as the option to make diagnostic plots.

    infile is an input_file object. See class InputFile

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
    See AGBTracks class.

    3. trilegal needs two files to link to Paola's formatted tracks
       a. track file goes here: trilegal_1.3/isotrack/tracce_[mix]_[set].dat
       b. cmd_input_file that links to the track file
          goes here trilegal_1.3/cmd_input_[mix]_[set].dat
    '''
    # set up file names and directories, cd to paola's tracks.
    AGB_file_setup(infile)

    # list of isofiles, zs, and ys to send to tracce file.
    isofiles, Zs, Ys = [], [], []
    imfr_data = np.array([])
    lifetime_data = np.array([])
    for metal_dir in infile.metal_dirs:
        metallicity, Y = metallicity_from_dir(metal_dir)
        print 'Z = %.4f' % metallicity
        if infile.diagnostic_dir0 is not None:
            diagnostic_dir = os.path.join(infile.diagnostic_dir0,
                                          infile.agb_mix,
                                          infile.set_name,
                                          metal_dir)  + '/'
            rsp.fileIO.ensure_dir(diagnostic_dir)
            # update infile class to place plots in this directory
            infile.diagnostic_dir = diagnostic_dir

        agb_tracks = rsp.fileIO.get_files(os.path.join(infile.working_dir,
                                                   metal_dir),
                                                   infile.track_identifier)
        agb_tracks.sort()

        iso_name_conv = '_'.join(('Z%.4f' % metallicity, infile.name_conv))
        isofile = os.path.join(infile.home, infile.isotrack_dir, iso_name_conv)

        if infile.over_write is False and os.path.isfile(isofile):
            print 'not over writing %s' % isofile
            out = None
        if infile.over_write is True:
            out = open(isofile, 'w')
            out.write('# age(yr) logL logTe m_act mcore c/o period ip')
            out.write(' Mdot(Msun/yr) logTe X Y CNO dlogTe/dlogL \n')

        isofile_rel_name = os.path.join('isotrack', infile.isotrack_dir,
                                        iso_name_conv)
        print 'found %i tracks' % len(agb_tracks)
        for agb_track in agb_tracks:
            # load track
            track = fileIO.get_numeric_data(agb_track)
            if track == -1:
                continue

            if track.bad_track is True:
                continue

            #if not track.mass in infile.masses:
            #    continue

            assert metallicity == track.metallicity, \
                'directory and track metallicity do not match'

            # make iso file for trilegal
            if out is not None:
                fileIO.make_iso_file(track, out)

            # save information for lifetime file.
            lifetime_datum = np.array([metallicity, track.mass, track.tauc,
                                       track.taum])

            lifetime_data = np.append(lifetime_data, lifetime_datum)

            # make diagnostic plots
            if infile.diagnostic_dir0 is not None and infile.diag_plots is True:
                assert metallicity_from_dir(infile.diagnostic_dir)[0] == \
                    track.metallicity, 'diag dir met wrong!'
                graphics.diag_plots(track, infile)

            # save information for imfr
            if infile.make_imfr is True:
                M_s = track.get_col('M_star')
                imfr_datum = np.array([M_s[0], M_s[-1], float(metallicity)])
                imfr_data = np.append(imfr_data, imfr_datum)

        if out is not None:
            out.close()
            print 'wrote', isofile

        fileIO.make_local_copy(isofile, dest=infile.make_copy)
        # keep information for tracce file
        isofiles.append(isofile_rel_name)
        Ys.append(Y)
        Zs.append(metallicity)
        #graphics.bigplots(agb_tracks, infile)

    # make file to link cmd_input to formatted agb tracks
    fileIO.make_met_file(infile.tracce_file, Zs, Ys, isofiles)

    fileIO.make_local_copy(infile.tracce_file, dest=infile.make_copy)

    # make cmd_input file
    cmd_in_kw = {'cmd_input_file': infile.cmd_input_file,
                 'file_tpagb': infile.tracce_file_rel,
                 'mass_loss': infile.mass_loss,
                 'file_isotrack': infile.file_isotrack}
    fileIO.write_cmd_input_file(**cmd_in_kw)

    if infile.make_imfr is True and infile.diagnostic_dir0 is not None:
        ifmr_file = os.path.join(infile.diagnostic_dir0, infile.agb_mix,
                                 infile.set_name,
                                 '_'.join(('ifmr', infile.name_conv)))
        ncols = 3
        nrows = imfr_data.size/ncols
        fileIO.savetxt(ifmr_file, imfr_data.reshape(nrows, ncols),
                       header='# M_i M_f Z \n')
        graphics.plot_ifmr(ifmr_file)

    if infile.diagnostic_dir0 is not None and infile.diag_plots is True:
        lifetime_file = os.path.join(infile.diagnostic_dir0, infile.agb_mix,
                                     infile.set_name,
                                     '_'.join(('tau_cm', infile.name_conv)))
        ncols = 4
        nrows = lifetime_data.size/ncols
        fileIO.savetxt(lifetime_file, lifetime_data.reshape(nrows, ncols),
                       header='# z mass tauc taum\n')
        graphics.plot_cluster_test(lifetime_file, infile)

    os.chdir(infile.home)
    return infile.cmd_input_file

def examine_1TP(agb_mix, set_name):
    outfile = 'examine1TP_%s_%s.dat' % (agb_mix, set_name)
    infile = 'cmd_input_%s_%s.dat' % (agb_mix, set_name)
    run_examine_1TP(infile, outfile=outfile)

def run_examine_1TP(cmd_inp_file, outfile=None):
    if outfile is None:
        outfile = '%s_examine1TP.dat' % cmd_inp_file
    here = os.getcwd()
    os.chdir(os.environ['TRILEGAL_ROOT'])
    os.system('./examine1TP.pl %s > %s' % (cmd_inp_file, outfile))
    lines = open(outfile, 'r').readlines()
    warns = [l for l in lines if 'Warning' in l]
    print 'Found %i warnings' % len(warns)
    print [l.strip() for l in warns]

    data = np.array([l for l in lines if l.startswith('Z')])
    nind = [i for i,l in enumerate(data) if not 'COLIBRI' in l]
    mass = np.array([m.strip().split()[3] for m in data], dtype=float)
    inds, = np.nonzero((mass>0.55) & (mass < 5.))
    missing = list(set(nind) & set(inds))
    print 'Missing TPAGB Tracks:'
    for i in missing:
        if len(data[i].split()) < len(data[i-1].split()):
            print data[i].strip()

    os.chdir(here)


if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
    except:
        print do_everything.__doc__
    pdb.set_trace()
    infile = rsp.fileIO.input_file(input_file, default_dict=fileIO.agb_input_defaults())

    agb_mix = infile.agb_mix
    set_name = infile.set_name
    track_set = '_'.join((agb_mix, set_name))


    # Paola's tracks -> trilegal + tests based only on tracks
    if infile.parse_tracks:
        do_everything(infile)
    if infile.examineAGB is True:
        examine_1TP(agb_mix, set_name)

    # Marco's scripts to run trilegal at age and z
    diagnostic_dir = infile.diagnostic_dir0
    if infile.trilegal_diagnostics and diagnostic_dir:
        tri_dir = os.path.join(diagnostic_dir, 'trilegal_files')
        sfh_dir = os.path.join(tri_dir, 'sfh')
        plt_dir = os.path.join(diagnostic_dir, agb_mix, set_name, track_set)
        trilegal_diagnostics.main(track_set, sfh_dir, tri_dir, plt_dir,
                                  over_write=infile.over_write, multi=False)
