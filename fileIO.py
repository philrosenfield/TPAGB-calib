import os
import numpy as np
import glob
import sys
import math_utils
import graphics

class AGBTracks(object):
    '''
    A fast way of loading files.

    examples:
    track = some agb track file name (string)
    AGB = get_numeric_data(track)

    all the data:
    AGB.data_array

    the file name:
    AGB.name

    a dictionary with 'key': column_number
    AGB.key_dict

    get the row of the AGB track by calling the number in the first column
    (ex: 7):
    AGB.get_row_bynum(7)

    get a specified column (use AGB.key_dict to know which one, or do head -1)
    (ex: logL)
    AGB.get_col('L_*')
    '''
    def __init__(self, data_array, col_keys, name):
        self.data_array = data_array
        self.key_dict = dict(zip(col_keys, range(len(col_keys))))
        self.name = name
        self.firstname = os.path.split(name)[1]
        self.mass = float(self.firstname.split('_')[1])
        self.metallicity = float(self.firstname.split('_')[2].replace('Z', ''))
        # initialize: it's a well formatted track with more than one pulse
        self.bad_track = False
        # if only one thermal pulse, stop the press.
        self.check_ntp()
        self.get_TP_inds()
        if not self.bad_track:
            # force the beginning phase to not look like it's quiescent
            self.fix_phi()
            # load quiescent tracks
            self.get_quiescent_inds()
            # load indices of m and c stars
            self.m_cstars()
            # calculate the lifetimes of m and c stars
            self.tauc_m()
            # add points to low mass quiescent track for better interpolation
            self.addpt = []
            if len(self.Qs) <= 9 and self.mass < 3.:
                self.add_points_to_q_track()
            # find dl/dt of track
            self.find_dldt()
        else:
            print 'bad track:', name

    def find_dldt(self, order=1):
        '''
        Finds dL/dt of track object by a poly fit of order = 1 (default)
        '''
        TPs = self.TPs
        qs = list(self.Qs)
        status = self.get_col('status')
        logl = self.get_col('L_star')
        logt = self.get_col('T_star')
        phi = self.get_col('PHI_TP')
        # if a low mass interpolation point was added it will get
        # the same slope as the rest of the thermal pulse.
        #dl/dT seems somewhat linear for 0.2 < phi < 0.4 ...
        lin_rise, = np.nonzero((status == 7) & (phi < 0.4) & (phi > 0.2))
        rising = [list(set(TP) & set(lin_rise)) for TP in TPs]
        fits = [np.polyfit(logt[r], logl[r], order) for r in rising if len(r) > 0]
        slopes = np.array([])
        # first line slope
        slopes = np.append(slopes, (logl[2] - logl[0]) / (logt[2] - logt[0]))
        # poly fitted slopes
        slopes = np.append(slopes, [fits[i][0] for i in range(len(fits))])

        # pop in an additional copy of the slope if an interpolation point 
        # was added.
        if len(self.addpt) > 0:
            tps_of_addpt = np.array([i for i in range(len(TPs))
                                    if list(set(self.addpt) & set(TPs[i])) > 0])
            slopes = np.insert(slopes, tps_of_addpt, slopes[tps_of_addpt])

        self.Qs = np.insert(self.Qs, 0, 0)
        self.rising = rising
        self.slopes = slopes
        self.fits = fits

    def add_points_to_q_track(self):
        '''
        when to add an extra point for low masses
        if logt[qs+1] is hotter than logt[qs]
        and there is a point inbetween logt[qs] and logt[qs+1] that is cooler
        than logt[qs] add the coolest point.
        '''
        addpt = self.addpt
        qs = list(self.Qs)
        logt = self.get_col('T_star')
        tstep = self.get_col('step')
        status = self.get_col('status')
        Tqs = logt[qs]
        # need to use some unique array, not logt, since logt could repeat,
        # index would find the first one, not necessarily the correct one.
        Sqs = tstep[qs] - 1.  # steps start at 1, not zero
        # takes the sign of the difference in logt(qs) 
        # if the sign of the difference is more than 0, we're going from cold to ho

        # finds where the logt goes from getting colder to hotter...
        ht, = np.nonzero(np.sign(np.diff(Tqs)) > 0)
        ht = np.append(ht, ht + 1)  # between the first and second
        Sqs_ht = Sqs[ht]
        # the indices between each hot point.
        t_mids = [map(int, tstep[int(Sqs_ht[i]): int(Sqs_ht[i + 1])])
                  for i in range(len(Sqs_ht) - 1)]
        Sqs_ht = Sqs_ht[: -1]
        for i in range(len(Sqs_ht) - 1):
            hot_inds = np.nonzero(logt[int(Sqs_ht[i])] > logt[t_mids[i]])[0]
            if len(hot_inds) > 0:
                # index of the min T of the hot index from above.
                addpt.append(list(logt).index(np.min(logt[[t_mids[i][hi]
                                                           for hi in hot_inds]])))

        if len(addpt) > 0:
            addpt = np.unique([a for a in addpt if status[a] == 7.])
        # hack: if there is more than one point, take the most evolved.
        if len(addpt) > 1:
            addpt = [np.max(addpt)]
        # update Qs with added pts.
        self.Qs = np.sort(np.concatenate((addpt, qs)))
        self.addpt = addpt

    def check_ntp(self):
        '''
        sets self.bad_track = True if only one thermal pulse.
        '''
        ntp = self.get_col('NTP')
        if ntp.size == 1:
            print 'no tracks!', self.name
            self.bad_track = True

    def fix_phi(self):
        '''
        The first line in the agb track is 1. This isn't a quiescent stage.
        '''
        self.data_array['PHI_TP'][0] = -1.

    def get_row(self, i):
        return self.data_array[i, :]

    def get_row_bynum(self, i):
        row = np.nonzero(self.data_array[:, self.key_dict['step']] == i)[0]
        return self.data_array[row, :]

    def get_col(self, key):
        return self.data_array[key]

    def m_cstars(self, mdot_cond=-5, logl_cond=3.3):
        '''
        adds mstar and cstar attribute of indices that are true for:
        mstar: co <=1 logl >= 3.3 mdot <= -5
        cstar: co >=1 mdot <= -5
        (by default) adjust mdot with mdot_cond and logl with logl_cond.
        '''
        data = self.data_array

        self.mstar, = np.nonzero((data['CO'] <= 1) &
                                 (data['L_star'] >= logl_cond) &
                                 (data['dMdt'] <= mdot_cond))
        self.cstar, = np.nonzero((data['CO'] >= 1) &
                                 (data['dMdt'] <= mdot_cond))

    def tauc_m(self):
        '''
        lifetimes of c and m stars
        '''
        try:
            tauc = np.sum(self.data_array['dt'][self.cstar]) / 1e6
        except IndexError:
            tauc = 0.
            print 'no tauc'
        try:
            taum = np.sum(self.data_array['dt'][self.mstar]) / 1e6
        except IndexError:
            taum = 0.
            print 'no taum'
        self.taum = taum
        self.tauc = tauc

    def get_TP_inds(self):
        '''
        find the thermal pulsations of each file
        '''
        if not self.bad_track:
            ntp = self.get_col('NTP')
            un = np.unique(ntp)
            if un.size == 1:
                print 'only one themal pulse.'
                self.TPs = un
            else:
                # this is the first step in each TP.
                iTPs = [list(ntp).index(u) for u in un]
                # The indices of each TP.
                TPs = [np.arange(iTPs[i], iTPs[i + 1]) for i in range(len(iTPs) - 1)]
                # don't forget the last one.
                TPs.append(np.arange(iTPs[i + 1], len(ntp)))
                self.TPs = TPs
        else:
            self.TPs = []
        if len(self.TPs) == 1:
            self.bad_track = True

    def get_quiescent_inds(self):
        '''
        The quiescent phase, Qs,  is the the max phase in each TP,
        i.e., closest to 1.
        '''
        phi = self.get_col('PHI_TP')
        self.Qs = np.unique([TP[np.argmax(phi[TP])] for TP in self.TPs])


def get_numeric_data(filename):
    '''
    made to read all of Paola's tracks. It takes away her "lg" meaning log.
    Returns an AGBTracks object. If there is a problem reading the data, all
    data are passed as zeros.
    '''
    f = open(filename, 'r')
    line = f.readline()
    f.close()
    col_keys = line.replace('#', '').replace('lg', '').replace('*', 'star')
    col_keys = col_keys.strip().split()
    try:
        data = np.genfromtxt(filename, missing_values='************',
                             names=col_keys)
    except ValueError:
        print 'problem with', filename
        data = np.zeros(len(col_keys))
    return AGBTracks(data, col_keys, filename)


def get_numeric_data_not_yet_working(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    col_keys = lines[0].replace('#', '').replace('lg', '').replace('*', 'star').replace('/','')
    col_keys = col_keys.strip().split()
    ncols = len(col_keys)
    nrows = len(lines) - 1
    dtype = [(c, '<f8') for c in col_keys]
    data = np.ndarray((nrows,), dtype=dtype)
    for i in range(len(lines)-1):
        if lines[i].startswith('#'):
            continue
        try:
            data[i] = np.array(map(float,lines[i].strip().split()))
        except ValueError:
            err = sys.exc_info()[1]
            tb = sys.exc_info()[-1]
            print ("%s %s"%(tb.tb_frame.f_code.co_filename, err))
    return AGBTracks(data, col_keys, filename)


def make_iso_file(track, isofile):
    '''
    this only writes the quiescent lines and the first line.
    format of this file is:
    t_min,          age in yr
    logl_min,       logL
    logte_min,      logTe
    mass_min,       actual mass along track
    mcore_min,      core mass
    co_min,         C/O ratio
    per_min,        period in days
    ip_min,         1=first overtone, 0=fundamental mode
    mlr_min,        - mass loss rate in Msun/yr
    logtem_min,     keep equal to logTe
    x_min,          X
    y_min,          Y
    xcno_min        X_C+X_O+X_N
    slope           dTe/dL
    '''
    fmt = '%.4e %.4f %.4f %.5f %.5f %.4f %.4e %i %.4e %.4f %.6e %.6e %.6e %.4f \n'

    # cull agb track to quiescent, write out.
    rows = [q for q in track.Qs]  # cutting the final point for the file.
    rows[0] += 1

    keys = track.key_dict.keys()
    vals = track.key_dict.values()
    col_keys = np.array(keys)[np.argsort(vals)]

    cno = [key for key in col_keys if (key.startswith('C1') or
                                       key.startswith('N1') or
                                       key.startswith('O1'))]

    isofile.write(' %.4f %i # %s \n' % (track.mass, len(rows), track.firstname))
    # hack, is the line length killing trilegal?
    # isofile.write(' %.4f %i # test \n' % (track.mass, len(rows)))
    for r in rows:
        row = track.data_array[r]
        CNO = np.sum([row[c] for c in cno])
        mdot = 10 ** (row['dMdt'])
        if row['Pmod'] == 0:
            period = row['P0']
        else:
            period = row['P1']
        if r == rows[-1]:
            # adding nonsense slope for the final row.
            slope = 999
        else:
            try:
                slope = 1. / track.slopes[list(rows).index(r)]
            except:
                print 'fucked',track.firstname
                graphics.hrd_slopes(track)
        try:
            isofile.write(fmt % (row['ageyr'], row['L_star'], row['T_star'],
                                 row['M_star'], row['M_c'], row['CO'], period,
                                 row['Pmod'], mdot, row['T_star'], row['H'],
                                 row['Y'], CNO, slope))
        except IndexError:
            print list(rows).index(r)
            print len(rows), len(track.slopes)
            print 1. / slopes[list(rows).index(r)]
    return


class input_file(object):
    '''
    a class to replace too many kwargs from the input file.
    does two things:
    1. sets a default dictionary (see input_defaults) as attributes
    2. unpacks the dictionary from load_input as attributes.
    '''
    def __init__(self, filename, default_dict=None):
        if default_dict is not None:
            self.set_defaults(default_dict)
        self.in_dict = load_input(filename)
        self.unpack_dict()

    def set_defaults(self, in_def):
        self.unpack_dict(udict=in_def)

    def unpack_dict(self, udict=None):
        if udict is None:
            udict = self.in_dict
        [self.__setattr__(k, v) for k, v in udict.items()]


def input_defaults(profile=None):
    '''
    the input file should be formatted in this way (comments are optional):
    # Paola's tracks are here:
    agbtrack_dir    /Users/phil/research/TP-AGBcalib/AGBTracks
    
    # Paola's directory structure is agbtrack_dir/agb_mix/set_name
    # This agb mix name
    agb_mix         CAF09    
    # This set name
    set_name        S_AUG12
    # Where cmd_input files will go for trilegal
    trilegal_dir    cmd_inputfiles/
    # any additional cmd_input parameters?
    mass_loss       Reimers: 0.35
    # Paola's formatted tracks are here: /[mix]/[set]
    isotrack_dir    isotrack_agb/
    # File to link from cmd_input to Paola's formatted tracks:
    tracce_dir      isotrack_agb/
    # make a another copy of the above files
    # (for example in dropbox or in a trilegal dir)
    make_copy       /Users/phil/research/padova_apps/isotrack_agb/
    # make initial and final mass relation (and also lifetimes c and m)?
    # (Need more than one metallicity)
    make_imfr       True
    # Diagnostic plots base 
    # this will have directory structure: 
    # diagnostic_dir0/[agb_mix]_[metallicity]/[set]/
    diagnostic_dir0            /Users/phil/research/TP-AGBcalib/diagnostics/
    # Run Marco's Scripts
    trilegal_diagnostics      True
    # Make google sites tables
    google_table    True
    # Do observational tests
    IDs   DDO82, NGC2403-HALO-6, NGC2976-DEEP, NGC4163, NGC7793-HALO-6, UGC8508, UGCA292
    # where to put the observational tests
    galaxy_outdir   /Users/phil/research/TP-AGBcalib/SNAP/models
    # overwrite Paola's parsed tracks?
    over_write   True
    # only do these metallicities (if commented out, do all metallicities)
    metals_subset       0.001, 0.004, 0.0005, 0.006, 0.008
    '''
    if profile is None:
        keys = ['over_write',
                'agbtrack_dir', 
                'agb_mix',
                'set_name',
                'trilegal_dir',
                'isotrack_dir',
                'tracce_dir',
                'diagnostic_dir0',
                'make_imfr',
                'make_copy',
                'mass_loss',
                'trilegal_diagnostics',
                'IDs',
                'galaxy_outdir',
                'image_location',
                'metals_subset']
    else:
        print 'only default profile set...'

    in_def = {}
    for k in keys:
        in_def[k] = None
    in_def['over_write'] = False
    return in_def


def write_cmd_input_file(**kwargs):
    '''
    make a TRILEGAL cmd_input file based on default.

    Send each parameter that is different than default by:
    kwargs = { 'kind_tpagb': 4, 'file_tpagb': 'isotrack/tracce_CAF09_S0.dat'}
    cmd_input_file = write_cmd_input_file(**kwargs)

    To make the default file:
    cmd_input_file = write_cmd_input_file()

    if you don't specify cmd_input_file, output goes to cmd_input_TEMP.dat
    '''
    kind_tracks = kwargs.get('kind_tracks', 2)
    file_isotrack = kwargs.get('file_isotrack', 'isotrack/parsec/CAF09.dat')
    file_lowzams = kwargs.get('file_lowzams', 'isotrack/bassazams_fasulla.dat')
    kind_tpagb = kwargs.get('kind_tpagb', 4)
    file_tpagb = kwargs.get('file_tpagb')
    if not file_tpagb:
        file_tpagb = 'isotrack/isotrack_agb/tracce_CAF09_AFEP02_I1_S1.dat'

    kind_postagb = kwargs.get('kind_postagb', 0)
    file_postagb = kwargs.get('file_postagb', 'isotrack/final/pne_wd_test.dat')
    mass_loss = kwargs.get('mass_loss')
    if mass_loss:
        kind_rgbmloss = 1
        law_mass_loss, = mass_loss.keys()
        efficiency_mass_loss, = mass_loss.values()
    # these are for using cmd2.2:
    kind_mag = kwargs.get('kind_mag', None)
    photsys = kwargs.get('photsys', 'wfpc2')
    file_mag = 'tab_mag_odfnew/tab_mag_%s.dat' % photsys
    kind_imf = kwargs.get('kind_imf', None)
    file_imf = kwargs.get('file_imf', 'tab_imf/imf_chabrier_lognormal.dat')

    # if not using cmd2.2:
    if kind_imf is None:
        kind_imfr = kwargs.get('kind_imfr', 0)
        file_imfr = kwargs.get('file_imfr', 'tab_ifmr/weidemann.dat')

    track_comments = '# kind_tracks, file_isotrack, file_lowzams'
    tpagb_comments = '# kind_tpagb, file_tpagb'
    pagb_comments = '# kind_postagb, file_postagb DA VERIFICARE file_postagb'
    mag_comments = '# kind_mag, file_mag'
    imf_comments = '# kind_imf, file_imf'
    imfr_comments = '# ifmr_kind, file with ifmr'
    mass_loss_comments = '# RGB mass loss: kind_rgbmloss, law, and efficiency'
    footer = (
        '################################explanation######################',
        'kind_tracks: 1= normal file',
        'file_isotrack: tracks for low+int mass',
        'file_lowzams: tracks for low-ZAMS',
        'kind_tpagb:',
        ' 0 = none',
        ' 1 = Girardi et al., synthetic on the flight, no dredge up',
        ' 2 = Marigo & Girardi 2001, from file, includes mcore and C/O',
        ' 3 = Marigo & Girardi 2007, from file, includes per, mode and mloss',
        ' 4 = Marigo et al. 2011, from file, includes slope',
        'file_tpagb: tracks for TP-AGB',
        'kind_postagb:',
        ' 0 = none',
        ' 1 = from file',
        'file_postagb: PN+WD tracks',
        'kind_ifmr:',
        ' 0 = default',
        ' 1 = from file\n')
    cmd_input_file = kwargs.get('cmd_input_file', 'cmd_input_TEMP.dat')
    fh = open(cmd_input_file, 'w')
    formatter = ' %i %s %s \n'
    fh.write(' %i %s %s %s \n' % (kind_tracks, file_isotrack, file_lowzams,
                                  track_comments))
    fh.write(formatter % (kind_tpagb, file_tpagb, tpagb_comments))
    fh.write(formatter % (kind_postagb, file_postagb, pagb_comments))
    if kind_mag is not None:
        fh.write(formatter % (kind_mag, file_mag, mag_comments))
    if kind_imf is None:
        fh.write(formatter % (kind_imfr, file_imfr, imfr_comments))
    else:
        fh.write(formatter % (kind_imf, file_imf, imf_comments))
    if mass_loss:
        fh.write(' %i %s %.3f \n' % (kind_rgbmloss, law_mass_loss,
                                     efficiency_mass_loss))
    fh.write('\n'.join(footer))
    fh.close()
    print 'wrote', cmd_input_file
    return cmd_input_file


def make_met_file(tracce, Zs, Ys, isofiles):
    with open(tracce, 'w') as t:
        t.write(' %i\n' % len(isofiles))
        [t.write(' %.4f\t%.3f\t%s\n' % (Zs[i], Ys[i], isofiles[i]))
         for i in np.argsort(Zs)]
    print 'wrote', tracce
    return
'''
def load_input(filename):
    in_def = input_defaults()
    infile = input_file(filename, default_dict=in_def)
    return infile
'''

def get_files(src, search_string):
    '''
    returns a list of files, similar to ls src/search_string
    '''
    if not src.endswith('/'):
        src += '/'
    try:
        files = glob.glob1(src, search_string)
    except IndexError:
        print 'Can''t find %s in %s' % (search_string, src)
        sys.exit(2)
    files = [os.path.join(src, f) for f in files]
    [ensure_file(f) for f in files] 
    return files


def load_input(filename):
    '''
    reads an input file into a dictionary.
    file must have key first then value(s)
    Will make 'True' into a boolean True
    Will understand if a value is a float, string, or list, etc.
    Ignores all lines that start with #, but not with # on the same line as
    key, value.
    '''
    try:
        literal_eval
    except NameError:
        from ast import literal_eval

    d = {}
    with open(filename) as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            if len(line.strip()) == 0:
                continue
            key, val = line.strip().partition(' ')[0::2]
            d[key] = math_utils.is_numeric(val.replace(' ', ''))
    # do we have a list?
    for key in d.keys():
        # float
        if type(d[key]) == float:
            continue
        # list:
        temp = d[key].split(',')
        if len(temp) > 1:
            try:
                d[key] = map(float, temp)
            except:
                d[key] = temp
        # dict:
        elif len(d[key].split(':')) > 1:
            temp1 = d[key].split(':')
            d[key] = {math_utils.is_numeric(temp1[0]): math_utils.is_numeric(temp1[1])}
        else:
            val = temp[0]
            # boolean
            true = val.upper().startswith('TR')
            false = val.upper().startswith('FA')
            if true or false:
                val = literal_eval(val)
            # string
            d[key] = val
    return d


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.isdir(d):
        os.makedirs(d)
        print 'made dirs: ', d


def savetxt(filename, data, fmt='%.4f', header=None):
    '''
    np.savetxt wrapper that adds header. Some versions of savetxt
    already allow this...
    '''
    with open(filename, 'w') as f:
        if header is not None:
            f.write(header)
        np.savetxt(f, data, fmt=fmt)
    print 'wrote', filename

def ensure_file(f, mad=True):
    '''
    input 
    f (string): if f is not a file will print "no file"
    optional
    mad (bool)[True]: if mad is True, will exit program.
    '''
    if not os.path.isfile(f):
        print 'there is no file', f
        if mad:
            sys.exit()

def make_local_copy(file, dest=os.environ['ISOTRACK'][:-1] + '_agb/'):
    if dest is not None:
        os.system('cp %s %s' % (file, dest))
    return
