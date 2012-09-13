import os
import numpy as np
import glob
import sys
import math_utils


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
        self.mass = float(os.path.split(name)[1].split('_')[1])
        firstname = os.path.split(name)[1]
        self.metallicity = float(firstname.split('_')[2].replace('Z', ''))

    def get_row(self, i):
        return self.data_array[i, :]

    def get_row_bynum(self, i):
        row = np.nonzero(self.data_array[:, self.key_dict['step']] == i)[0]
        return self.data_array[row, :]

    def get_col(self, key):
        return self.data_array[key]


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

def make_iso_file(track, Qs, slopes, isofile):
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
    fmt = (' %.4e %.4f %.4f %.5f %.5f %.4f %.4e %i %.4e %.4f %.6e %.6e %.6e %.4f\n')

    # cull agb track to quiescent, write out.
    rows = [q - 1 for q in Qs]  # cutting the final point for the file.
    rows[0] += 2

    keys = track.key_dict.keys()
    vals = track.key_dict.values()
    col_keys = np.array(keys)[np.argsort(vals)]

    cno = [key for key in col_keys if (key.startswith('C1') or
                                       key.startswith('N1') or
                                       key.startswith('O1'))]

    isofile.write(' %.4f %i # %s\n' % (track.mass, len(rows),
                                       os.path.split(track.name)[1]))
    for r in rows:
        data_row = track.data_array[r]
        CNO = np.sum([data_row[c] for c in cno])
        mdot = 10 ** (data_row['dMdt'])
        if data_row['Pmod'] == 0:
            period = data_row['P0']
        else:
            period = data_row['P1']
        if r == rows[-1]:
            # adding nonsense slope for the final row.
            slope = 999999
        else:
            slope = 1. / slopes[list(rows).index(r)]
        try:
            isofile.write(fmt % (data_row['ageyr'], data_row['L_star'],
                                 data_row['T_star'], data_row['M_star'],
                                 data_row['M_c'], data_row['CO'], period,
                                 data_row['Pmod'], mdot, data_row['T_star'],
                                 data_row['H'], data_row['Y'], CNO, slope))
        except IndexError:
            print list(rows).index(r)
            print len(rows), len(slopes)
            print 1. / slopes[list(rows).index(r)]
    return


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
    return cmd_input_file


def make_met_file(tracce, Zs, Ys, isofiles):
    t = open(tracce, 'w')
    t.write(' %i\n' % len(isofiles))
    [t.write(' %.4f\t%.3f\t%s\n' % (Zs[i], Ys[i], isofiles[i]))
     for i in np.argsort(Zs)]
    t.close()
    return


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.isdir(d):
        os.makedirs(d)
        print 'made dirs: ', d


def get_files(src, search_string):
    '''
    simple search, returns a list
    '''
    import glob
    try:
        files = glob.glob1(src, search_string)
    except IndexError:
        print 'Can''t find', search_string, 'in', src
        sys.exit(2)
    return [os.path.join(src, f) for f in files]


def make_local_copy(file, dest=os.environ['ISOTRACK'][:-1] + '_agb/'):
    if dest is not None:
        os.system('cp %s %s' % (file, dest))
    return
