import os
import sys
import glob
try:
    import numpy as np
except ImportError:
    print 'You need numpy to run this script'
    print 'sudo apt-get install python-numpy'
    sys.exit()

if np.__version__ < '1.6':
    print 'Warning -- you should update numpy, you might get errors. E.g:'
    print 'sudo apt-get install python-numpy=1:1.7.1-1ubuntu1'


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
    lines = f.readlines()
    f.close()
    line = lines[0]
    if len(lines) == 1:
        print 'only one line in %s' % filename
        return -1
    col_keys = line.replace('#', '').replace('lg', '').replace('*', 'star')
    col_keys = col_keys.strip().split()
    try:
        data = np.genfromtxt(filename, missing_values='************',
                             names=col_keys)
    except ValueError:
        print 'problem with', filename
        data = np.zeros(len(col_keys))
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
    # I don't know why this is here.
    #rows[0] += 1

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
                print 'SLOPE ERROR', track.firstname
        try:
            isofile.write(fmt % (row['ageyr'], row['L_star'], row['T_star'],
                                 row['M_star'], row['M_c'], row['CO'], period,
                                 row['Pmod'], mdot, row['T_star'], row['H'],
                                 row['Y'], CNO, slope))
        except IndexError:
            print list(rows).index(r)
            print len(rows), len(track.slopes)
            print 1. / slope[list(rows).index(r)]
    return


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
    return files


def metallicity_from_dir(met):
    ''' take Z and Y values from string'''
    if met.endswith('/'):
        met = met[:-1]

    if len(os.path.split(met)) > 0:
        met = os.path.split(met)[1]

    z = float(met.split('_')[1].replace('Z', ''))
    y = float(met.split('_')[-1].replace('Y', ''))

    return z, y


def parse_tracks(metal_dir, track_identifier='agb_*Z*.dat',
                  agb_mix='CAF09', set_name='AGBSET'):
    '''
    parse colibri agb tracks.

    parameters:
    metal_dir string
        the prefix of the directory with agb tracks, e.g., S12_Z0.0001_Y0.249
    track_identifier string
        a unix like search string for the agb track, do not set to agbq.
    agb_mix string
        the agb mix for the file name e.g, CAF09
    set_name string
        this agb model for the file name e.g., S_NOV13

    returns:
    writes parsed files to directory ./isotrack_agb

    usage:
    e.g. parse the tracks of CAF09/S_NOV13/S12_Z0.001_Y0.250:
    navigate to CAF09/S_NOV13/

    $ python parse_agb_tracks.py CAF09 S_NOV13 S12_Z0.001_Y0.250
    '''
    # the current directory
    working_dir = os.getcwd()

    # find the metallicity
    metallicity, Y = metallicity_from_dir(metal_dir)
    print 'Z = %.4f' % metallicity
    # load the list of agb tracks and sort
    agb_tracks = get_files(os.path.join(working_dir, metal_dir),
                           track_identifier)
    agb_tracks.sort()
    print 'found %i tracks' % len(agb_tracks)

    # set up and begin writting the parsed file
    name_conv = '%s.dat' % '_'.join(('Z%.4f' % metallicity, agb_mix, set_name))

    # make sure we can write to the isotrack_agb directory
    ensure_dir(os.path.join(working_dir, 'isotrack_agb/'))

    isofile = os.path.join(working_dir, 'isotrack_agb', name_conv)
    isofile_header = '# age(yr) logL logTe m_act mcore c/o period ip'
    isofile_header += ' Mdot(Msun/yr) logTe X Y CNO dlogTe/dlogL \n'
    out = open(isofile, 'w')
    out.write(isofile_header)

    for agb_track in agb_tracks:
        # load track
        track = get_numeric_data(agb_track)

        # skip bad tracks
        if track == -1:
            continue
        if track.bad_track is True:
            continue

        # parse the track
        make_iso_file(track, out)

    out.close()
    print 'wrote %s' % isofile
    return


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.isdir(d):
        os.makedirs(d)
        print 'made dirs: ', d

if __name__ == "__main__":
    try:
        agb_mix, set_name, metal_dir = sys.argv[1:]
    except:
        print parse_tracks.__doc__
    parse_tracks(metal_dir, agb_mix=agb_mix, set_name=set_name)

