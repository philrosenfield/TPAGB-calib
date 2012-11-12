from optparse import OptionParser
import operator
import cmdUtils
import brewer2mpl
import ResolvedStellarPops as rsp
from TrilegalUtils import get_stage_label
import LFUtils
import mk_sims
from multiprocessing import Pool
import itertools
import time
from TPAGBparams import *
import numpy as np
#import pdb; pdb.set_trace()
import matplotlib.pyplot as plt
import GoogleSitesTable as gst
import MatchUtils
import logging
logger = logging.getLogger()


def edit_sfh(ID, new_sfr_dir, orig_sfr_dir=None, **changes):
    '''
    add or subtract a value from each bin in sfr file.

    possible changes:
    age
    sfr
    z
    feh (dex)
    zdisp - - not yet.

    example:
    IDs = all_IDs()
    new_sfr_dir = os.path.join(data_src, 'sfh - 0.3dex')
    new_sfrs = [edit_sfh(ID, new_sfr_dir, **{'FeH': - 0.3}) for ID in IDs]
    '''
    z_sun = 0.01524
    sfr_file = mk_sims.get_sfrFILE(ID, sfr_dir=orig_sfr_dir)
    orig_sfr = np.genfromtxt(sfr_file, names=['age', 'sfr', 'z'])
    new_sfr = orig_sfr.copy()
    for key in changes.keys():
        if key in new_sfr.dtype.names:
            new_sfr[key] += changes[key]
        elif key.lower() == 'feh':
            z = new_sfr['z']
            feh = np.log10(z / z_sun)
            new_feh = feh + changes[key]
            new_sfr['z'] = z_sun * 10 ** new_feh
        else:
            print '%s not found in %s' % (key, sfr_file)
            return -1

    rsp.fileIO.ensure_dir(new_sfr_dir)
    new_sfr_file = os.path.join(new_sfr_dir, os.path.split(sfr_file)[1])
    np.savetxt(new_sfr_file, new_sfr, fmt='%2.6f')
    return new_sfr_file


def tag_cmds(IDs):
    '''
    only need to do this once, hopefully. Just went though and added integer
    tags to each fits table so that it would be like trilegal output with the
    - l option.
    '''
    seqs = ['MS', 'RHeB', 'RGB']
    if type(IDs) == str:
        search_string = '*' + '*'.join((IDs, 'IR', '.fits'))
        fits = rsp.fileIO.get_files(fits_src, search_string)
    else:
        search_strings = ['*' + '*'.join((ID, 'IR', '.fits')) for ID in IDs]
        fits = np.squeeze([rsp.fileIO.get_files(fits_src, ss)
                           for ss in search_strings])
    rsp.annotate_cmd.define_color_mag_region(fits, seqs)
    return


def run_all(IDs, models):
    pool = Pool()

    res = []
    for ID, model in itertools.product(IDs, models):
        res.append(pool.apply_async(main, (ID, model),))

    for r in res:
        r.get()

    return


def read_tagged_phot(tagged_file):
    '''
    reads a file created by rsp.annotate_cmd.define_color_mag_region. 
    ascii with 7 columns.
    '''
    if type(tagged_file) == str:
        #print 'reading %s' % tagged_file
        cols = ['ra', 'dec', 'mag1', 'mag2', 'mag1err', 'mag2err', 'stage']
        fits = np.genfromtxt(tagged_file, names=cols)
    else:
        fits = tagged_file
    return fits


def get_trgb_fitsname(ID, band):
    if band == 'opt':
        fits = rsp.fileIO.get_files(fits_src,
                                  '*' + '*'.join((ID, 'trim', '.fits')))[0]
        trgb = get_tab5_trgb_Av_dmod(ID)[0]
    elif band == 'ir':
        # read data
        fits, = rsp.fileIO.get_files(fits_src,
                                   '*' + '*'.join((ID, 'IR', '.fits')))
        trgb = LFUtils.get_trgb_ir_nAGB(ID)[0]
    else:
        print 'choose opt or ir'
        sys.exit()
    return trgb, fits


class galaxies(object):
    '''
    THIS IS FROM BR RATIO CODE. I'M NOT SURE IT'S GENERAL ENOUGH TO PUT INTO
    ~ / research / python SO I COPIED IT...

    I made summary...
    wrapper for lists of galaxy objects, each method returns lists, unless they
    are setting attributes.
    '''
    def __init__(self, galaxy_objects):
        self.galaxies = galaxy_objects
        galaxies.summary(self, 'target', 'filter1', 'filter2')
        # maybe I can avoid simgalaxies with self.__class__?
        #self.zs = np.unique([np.round(g.z, 4) for g in self.galaxies])
        #self.photsyss =  np.unique(g.photsys for g in self.galaxies)

    def summary(self, *attrs):
        '''
        similar to squish, but has np.unique.
        '''
        for attr in attrs:
            new_list = [g.__getattribute__(attr) for g in self.galaxies]
            new_val = np.unique(new_list)
            self.__setattr__('%ss' % attr, new_val)

    def all_stages(self, *stages):
        '''
        adds the indices of any stage as attributes to galaxy.
        If the stage isn't found, - 1 is returned.
        '''
        [g.all_stages(*stages) for g in self.galaxies]
        return

    def squish(self, *attrs):
        '''
        concatenates an attribute or many attributes and adds them to galaxies
        instance.
        No slicing, so not sure how it will be useful besides Color Mag2.
        '''
        for attr in attrs:
            new_list = [g.__getattribute__(attr) for g in self.galaxies]
            new_val = np.concatenate(new_list)
            self.__setattr__('%ss' % attr, new_val)

    def finite_key(self, key):
        return [g for g in self.galaxies if np.isfinite(g.__dict__[key])]

    def select_on_key(self, key, val):
        ''' ex filter2 == F814W works great with strings or exact g.key == val.
        rounds z to four places, no error handling.
        '''
        key = key.lower()
        if key == 'z':
            gs = [g for g in self.galaxies if
                  np.round(g.__dict__[key], 4) == val]
        else:
            gs = [g for g in self.galaxies if g.__dict__[key] == val]
        return gs

    def group_by_z(self):
        zsf = self.zs[np.isfinite(self.zs)]
        zsn = self.zs[np.isnan(self.zs)]
        d = {}
        for z in zsf:
            key = 'Z%.4f' % z
            d[key] = galaxies.select_on_key(self, 'z', z)

        d['no z'] = [g for g in gals if np.isnan(g.z)]
        return d

    def intersection(self, **kwargs):
        '''
        ex kwargs = {'filter2': 'F814W', 'filter1': 'F555W'}
        will return a list of galaxy objects that match all kwarg values.
        '''
        gs_tmp = self.galaxies
        gs = [galaxies.select_on_key(self, k, v) for k, v in kwargs.items()]
        for i in range(len(gs)):
            gs_tmp = list(set(gs_tmp) & set(gs[i]))
        return gs_tmp

    def __str__(self):
        for g in self.galaxies:
            print g.__str__()
        return ''


class galaxy(object):
    def __init__(self, ID, band):
        self.name = get_trgb_fitsname(ID, band)[1]
        self.data = read_tagged_phot(rsp.fileIO.replace_ext(self.name, '.dat'))
        self.base = os.path.split(self.name)[0]
        self.target = cmdUtils.get_fileinfo(self.name)[0]
        self.filter1 = cmdUtils.get_fileinfo(self.name)[1]
        self.filter2 = cmdUtils.get_fileinfo(self.name)[2]
        self.mag1 = self.data['mag1']
        self.mag2 = self.data['mag2']
        self.color = self.data['mag1'] - self.data['mag2']
        self.stage = self.data['stage']
        self.trgb = get_trgb_fitsname(ID, band)[0]
        self.z = LFUtils.get_key_fromtable(ID, 'Z')

    def stage_inds(self, stage_name):
        return get_stage_inds(self.data, stage_name)

    def delete_data(self):
        data_names = ['data', 'mag1', 'mag2', 'color', 'stage']
        [self.__delattr__(data_name) for data_name in data_names]

    def all_stages(self, *stages):
        '''
        adds the indices of some stage as an attribute.
        '''
        for stage in stages:
            i = stage_inds(self.stage, stage)
            self.__setattr__('i%s' % stage.lower(), i)
        return


class simgalaxy(object):
    def __init__(self, trilegal_out, filter1, filter2):
        self.base, self.name = os.path.split(trilegal_out)
        self.data = rsp.fileIO.read_table(trilegal_out)
        self.filter1 = filter1
        self.filter2 = filter2
        self.target = self.name.split('_')[2]
        self.mag1 = self.data.get_col(self.filter1)
        self.mag2 = self.data.get_col(self.filter2)
        self.stage = self.data.get_col('stage')

        simgalaxy.load_ast_corrections(self)

        data_to_slice = ['mag1', 'mag2', 'stage', 'ast_mag1', 'ast_mag2']
        slice_inds = self.rec
        simgalaxy.slice_data(self, data_to_slice, slice_inds)

        self.ast_color = self.ast_mag1 - self.ast_mag2
        self.color = self.mag1 - self.mag2
        simgalaxy.load_ic_mstar(self)

    def get_fits(self):
        match_out_dir = os.path.join(os.path.split(self.base)[0], 'match',
                                     'output')
        fit_file_name = '%s_%s_%s.fit' % (self.ID, self.mix, self.model_name)
        try:
            fit_file, = fileIO.get_files(match_out_dir, fit_file_name)
            self.chi2, self.fit = MatchUtils.get_fit(fit_file)
        except ValueError:
            print 'no match output for %s.' % fit_file_name
        return

    def load_ast_corrections(self):
        diff1 = self.data.get_col('diff_' + self.filter1)
        diff2 = self.data.get_col('diff_' + self.filter2)
        recovered1, = np.nonzero(abs(diff1) < 90.)
        recovered2, = np.nonzero(abs(diff2) < 90.)
        self.rec = list(set(recovered1) & set(recovered2))
        self.ast_mag1 = self.mag1 + diff1
        self.ast_mag2 = self.mag2 + diff2

    def slice_data(self, data_to_slice, slice_inds):
        '''
        slice already set attributes by some index list.
        '''
        [self.__setattr__(d, self.__dict__[d][slice_inds])
         for d in data_to_slice]

    def mix_modelname(self, model):
        self.mix, self.model_name = get_mix_modelname(model)

    def delete_data(self):
        '''
        for wrapper functions, I don't want gigs of data stored when they
        are no longer needed.
        '''
        data_names = ['data', 'mag1', 'mag2', 'color', 'stage', 'ast_mag1',
                      'ast_mag2', 'ast_color', 'rec']
        [self.__delattr__(data_name) for data_name in data_names]

    def stage_inds(self, stage_name):
        return np.nonzero(self.stage == get_stage_label(stage_name))[0]

    def load_ic_mstar(self):
        co = self.data.get_col('C/O')[self.rec]
        lage = self.data.get_col('logAge')[self.rec]
        mdot = self.data.get_col('logML')[self.rec]
        logl = self.data.get_col('logL')[self.rec]

        self.imstar, = np.nonzero((co <= 1) &
                                  (logl >= 3.3) &
                                  (mdot <= - 5) &
                                  (self.stage == get_stage_label('TPAGB')))

        self.icstar, = np.nonzero((co >= 1) &
                                  (mdot <= - 5) &
                                  (self.stage == get_stage_label('TPAGB')))

    def all_stages(self, *stages):
        '''
        adds the indices of some stage as an attribute.
        '''
        for stage in stages:
            i = stage_inds(self.stage, stage)
            self.__setattr__('i%s' % stage.lower(), i)
        return


def get_mix_modelname(model):
    mix = model.replace('cmd_input','').split('.')[0].split('_')[0]
    model_name = '_'.join(model.split('.')[0].split('_')[1:])
    return mix, model_name


def match_metallicities(IDs):
    sfhs = [mk_sims.get_sfrFILE(ID) for ID in IDs]
    zs = {}
    for ID, sfh in zip(IDs, sfhs):
        age, sfr, z = np.loadtxt(sfh, unpack=True)
        zs[ID] = {'z': z, 'avez': np.mean(z)}
    return zs


def compare_metallicities():
    IDs = allIDs()
    zs = match_metallicities(IDs)
    for ID in IDs:
        zs[ID]['zmeas'] = LFUtils.get_key_fromtable(ID, 'Z')

    for id, zdict in sorted(zs.items(), key=lambda(k, v): (v['avez'], k)):
        print '%s %.4f %.4f' % (id, zdict['avez'], zdict['zmeas'])


def load_galaxy_tagged(ID, band):
    trgb, fitsname = get_trgb_fitsname(ID, band)
    tagged_fits = read_tagged_phot(rsp.fileIO.replace_ext(fitsname, '.dat'))

    target, filt1, filt2 = cmdUtils.get_fileinfo(fitsname)
    mag1 = tagged_fits['mag1']
    mag2 = tagged_fits['mag2']
    stage = tagged_fits['stage']

    color = mag1 - mag2

    return tagged_fits


def get_stage_inds(fits, stage_name):
    inds, = np.nonzero(fits['stage'] == get_stage_label(stage_name))
    return inds


def make_normalized_simulation(ID, model, **kwargs):
    '''
    Wrapper for running trilegal and LFUtils.calc_LF
    Will continue to run trilegal until galaxy is high enough mass to have
    proper normalization.

    many attributes of galaxy and simgalaxy instance are set in
    LFUtils.calc_LF.

    inputs:
    ID of target: string
    cmd_input model to run trilegal: string

    optional params:
    over_write: bool [False]
        force running of trilegal (will run if no output file exists)
    maglims: list, tuple or string. ['trgb']
        either the dim, bright maglims or 'trgb'
    offsets: list, tuple: (dim, bright) default: [(1.5, 0.5)]
        must exist if maglims='trgb'

        maglims will be trgb + dim, trgb + bright
    object_mass: float [5e6]
        initial mass to run trilegal.

    returns:
        gal instance
        sgal instance
        ks test p_value
        maglims
    '''
    # To be compatible with Jason, might work with 'ir' or 'opt' should only
    # call different file finding routines.
    band = kwargs.get('band', 'ir')
    # load kwargs
    over_write = kwargs.get('over_write', False)
    object_mass = kwargs.get('object_mass', 5e6)
    maglims = kwargs.get('maglims', 'trgb')
    offsets = kwargs.get('offsets', (1.5, 0.5))
    mk_sims_args = kwargs.get('mk_sims_args', {})

    # data galaxy instance
    gal = galaxy(ID, band)
    # set maglim if not already set.
    if maglims == 'trgb':
        maglims = (gal.trgb + offsets[0], gal.trgb + offsets[1])

    # initializations
    go = 0
    normalization = 1e9
    # arbitrary normalization factor set to 0.75... 1.0 is upper limit,
    # didn't think it's necessary to make this a kwarg...
    while normalization > .75:
        # object_mass increase factor is hard coded, could be kwarg.
        if go != 0:
            object_mass = object_mass * 5.
            # reset over_write for while loop
            over_write = 1
        go += 1

        # run trilegal
        if over_write:
            print 'Trying %s %s, Mass = %g, Attempt %i' % (ID, model,
                                                           object_mass, go)
        mk_sims_args['object_mass'] = object_mass
        mk_sims_args['over_write'] = over_write
        trilegal_out = mk_sims.mk_sims(ID, model, **mk_sims_args)

        # load simgalaxy
        sgal = simgalaxy(trilegal_out, gal.filter1, gal.filter2)
        sgal.ID = ID
        sgal.model = model
        sgal.mix_modelname(model)

        # calcuate the LFs
        p_value = LFUtils.calc_LF(gal, sgal, maglims)

        # update normalization
        normalization = sgal.normalization
        print 'normalization', normalization

    if over_write:
        sgal.object_mass = object_mass
        print '\n %s necessary object_mass = %g\n' % (gal.name, object_mass)

    return gal, sgal, p_value, maglims


def all_IDs():
    IDs = ["SCL-DE1",
           "NGC2403-HALO-6",
           "NGC7793-HALO-6",
           "UGC-04459",
           "UGC-4305-1",
           "UGC-4305-2",
           "UGC-5139",
           "DDO78",
           "DDO82",
           "KDG73",
           "KKH37",
           "NGC0300-WIDE1",
           "NGC404",
           "NGC2403-DEEP",
           "NGC2976-DEEP",
           "NGC3741",
           "NGC4163",
           "UGC8508",
           "UGCA292",
           "NGC3077-PHOENIX",
           "IC2574-SGS",
           "HS117",
           "DDO71"]
    return IDs


def main(IDs=None, models=None, **kwargs):
    '''
    big wrapper for galaxy tests

    kwargs:
    make_plots: LFUtils.plotLFIR and .diagnostic_cmd
    publish_plots: call GoogleSitesTable
    '''
    make_plots = kwargs.get('make_plots')
    publish_plots = kwargs.get('publish_plots')
    sfr_dir = kwargs.get('sfrdir')
    out_dir = kwargs.get('outdir')
    if IDs is None:
        print 'give me a galaxy!'
        return -1

    if models is None:
        print 'give me a model!'
        return -1

    outfile = kwargs.get('outfile', 'result_tab.dat')
    if os.path.isfile(outfile):
        print outfile, 'exists. appending.'
        out = open(outfile, 'a')
    else:
        out = open(outfile, 'w')
        out.write('# ID model p_value NRGB_data NAGB_data NRGB_model NAGB_model\n')

    for ID, model in itertools.product(IDs, models):
        mk_sims_args = {'sfr_dir': sfr_dir, 'outdir': out_dir}
        sim_kwargs = {'mk_sims_args': mk_sims_args}
        gal, sgal, p_value, maglims = make_normalized_simulation(ID, model,
                                                                 **sim_kwargs)
        write_spread_catalog(sgal, **mk_sims_args)

        if make_plots:
            fig_loc = os.path.join(plt_dir, sgal.mix, sgal.model_name)
            diag_loc = os.path.join(fig_loc, 'diag')
            rsp.fileIO.ensure_dir(diag_loc + '/')

            figname = os.path.join(diag_loc,
                                   '%s_%s_diag.png' % (ID, sgal.model_name))
            LFUtils.diagnostic_cmd(sgal, gal.trgb, figname=figname,
                                   inds=sgal.rel_ind)
            # should take fig_loc...
            LFUtils.plot_LFIR(gal, sgal, p_value, maglims)

        out.write('%s %s %.3f %i %i %i %i\n' % (sgal.ID,
                                                sgal.model,
                                                p_value,
                                                gal.irgb.size,
                                                gal.iagb.size,
                                                sgal.rel_rgb.size,
                                                sgal.rel_agb.size))

    out.close()

    for model in models:
        #compare_sims(IDs, model)
        if publish_plots:
            html_file = os.path.join(diag_loc, sgal.model_name + '_diag.html')
            gst.side_by_side(diag_loc, sgal.model_name, html_file)

            html_file = os.path.join(fig_loc, sgal.model_name + '.html')
            gst.one_col(fig_loc, sgal.model_name, html_file)

    return


def load_galaxies(IDs, model, **kwargs):
    sgals = []
    gals = []
    #IDs = all_IDs()
    out_dir = kwargs.get('out_dir', model_src)
    #for ID, model in itertools.product(IDs, models):
    for ID in IDs:
        spread_file = os.path.join(out_dir, 'spread',
                                   'spread_output_%s_model_%s' % (ID, model))
        gal = galaxy(ID, 'ir')
        sgal = simgalaxy(spread_file, gal.filter1, gal.filter2)
        sgal.ID = ID
        sgal.model = model
        sgal.mix_modelname(model)
        sgal.get_fits()
        maglims = (np.nan, np.nan)
        p_value = LFUtils.calc_LF(gal, sgal, maglims, normalize=False)
        sgal.file_ast = mk_sims.get_fakFILE(ID)
        sgals.append(sgal)
        gals.append(gal)

    # galaxies - ify
    SGals = galaxies(sgals)
    Gals = galaxies(gals)
    return Gals, SGals


def load_galaxies_by_z(IDs, model):
    Gals, SGals = load_galaxies(IDs, model)
    # select on the ones with measured z:
    galsz = Gals.finite_key('z')
    # sort by z.
    galsz_sort = sorted(galsz, key=lambda galaxy: galaxy.z)

    # use that to grab the sim gals
    sgalsz = []
    for g in galsz_sort:
        sgalsz += [sg for sg in SGals.galaxies if sg.target == g.target]

    zs = [g.z for g in galsz_sort]
    return zs, galsz_sort, sgalsz


def compare_sims(IDs, model, fig_name=None):
    zs, galsz_sort, sgalsz = load_galaxies_by_z(IDs, model)
    ax = plt.axes()
    for z, g, s in zip(zs, galsz_sort, sgalsz):
        ax.plot(z, float(s.rel_agb.size) / float(g.iagb.size),
                'o', ms=5, color='black')
        bright_rgb = LFUtils.brighter(s.ast_mag2, g.trgb, inds=s.irgb)
        ax.plot(z,
                float((s.rel_agb.size - bright_rgb.size)) / float(g.iagb.size),
                'o', ms=5, color='red')
        ax.annotate(g.target,
                    xy=(z, float((s.rel_agb.size - bright_rgb.size)) /
                        float(g.iagb.size) + 0.01), ha='center')
        ax.annotate(g.target,
                    xy=(z, float(s.rel_agb.size) / float(g.iagb.size) + 0.01),
                    ha='center')
    ax.set_ylabel(r'$N_{AGB, model} / N_{AGB, data}$', fontsize=20)
    ax.set_xlabel(r'$Z$', fontsize=20)
    ax.set_title(r'$\rm{%s}$' % s.model_name.replace('_', '\ '))
    rsp.fileIO.ensure_dir(os.path.join(plt_dir, s.mix))
    if not fig_name:
        fig_name = os.path.join(plt_dir, s.mix, s.model_name + '.png')
    plt.savefig(fig_name)
    plt.close()
    print 'wrote %s' % fig_name
    return


def chi2_plot(IDs, models):
    ax = plt.axes()
    cols = brewer2mpl.get_map('Dark2', 'qualitative', len(models)).mpl_colors
    l = 0.
    for i, model in enumerate(models):
        l += .004
        zs, galsz_sort, sgalsz = load_galaxies_by_z(IDs, model)
        for z, g, s in zip(zs, galsz_sort, sgalsz):
            ax.plot(z, s.chi2, 'o', ms=5, color=cols[i])
            ax.annotate(g.target, xy=(z + 0.0001, s.chi2 + 0.01), ha='left',
                        va='center', fontsize='8')
        ax.annotate('$%s$' % s.model_name.replace('_', '\ '), xy=(l, 3),
                    color=cols[i], fontsize=20)
    ax.set_ylabel(r'$\chi^2$', fontsize=20)
    ax.set_xlabel(r'$Z$', fontsize=20)
    fig_name = os.path.join(plt_dir, s.mix, 'chi2.png')
    ax.set_ylim(2, 18)
    plt.savefig(fig_name)
    plt.close()
    return


def write_spread_catalog(sgal, outfile=None, **kwargs):
    '''
    writes a slice of the trilegal output catalog.
    slice is on ast recovered stars that are randomly selected by
    normalization.
    i.e sgal.rec and sgal.rel_ind.
    returns string outfile name.
    '''
    out_dir = kwargs.get('outdir')
    if out_dir is None:
        out_dir = model_src
    if outfile is None:
        outfile = os.path.join(out_dir, 'spread', sgal.name.replace('ast',
                                                                    'spread'))
        rsp.fileIO.ensure_dir(outfile)
    if not os.path.isfile(outfile):
        out = open(outfile, 'w')
        # get the header
        sort_keys = sorted(sgal.data.key_dict.iteritems(),
                           key=operator.itemgetter(1))
        col_keys = [k[0] for k in sort_keys]

        # slice data on recovered inds and rel_inds.
        dslice = np.squeeze(sgal.data.data_array[[sgal.rec], :])[sgal.rel_ind]

        out.write('# {}\n'.format(' '.join(col_keys)))
        np.savetxt(out, dslice, fmt='%-7.6f')
        print 'wrote %s' % outfile
    else:
        print '%s exists. Write spread catalog will not overwrite.' % outfile

    return outfile


def setup_match(name_fmt, bg_name=None, phot_name=None, **kwargs):
    '''
    sets up the match directory structure
    '''
    if bg_name is not None and phot_name is not None:
        mdirs = {'msgs': '.msg',
                 'output': '.fit',
                 'plots': '',
                 'pars': '.par',
                 'bg': bg_name,
                 'phot': phot_name}
        match_dir = os.path.join(model_src, 'match')
    else:
        match_dir = kwargs.get('match_dir')
        kwargs.pop('match_dir')
        mdirs = kwargs

    mdirdict = {}
    for d, ext in mdirs.items():
        new_dir = os.path.join(match_dir, d + '/')
        rsp.fileIO.ensure_dir(new_dir)
        # if no extension, return the plot directory
        if len(ext) == 0:
            mdirdict[d] = new_dir
            continue
        # if extension less than 5 chars, use name_fmt.ext
        if len(ext) <= 5:
            mdirdict[d] = os.path.join(new_dir, '%s%s' % (name_fmt, ext))
        # if more than 5 characters, assume it's a file name.
        else:
            mdirdict[d] = os.path.join(new_dir, ext)

    return mdirdict


def match_tests(IDs, model):
    Gals, SGals = load_galaxies(IDs, model)
    for g, s in zip(Gals.galaxies, SGals.galaxies):
        name_fmt = '%s_%s_%s' % (s.ID, s.mix, s.model_name)
        bg_name = s.name.replace('spread', 'bg')
        phot_name = rsp.fileIO.replace_ext(os.path.split(g.name)[1], '.match')
        match_dict = setup_match(name_fmt, bg_name, phot_name)
        match_dict['plots'] = os.path.join(match_dict['plots'], s.mix,
                                           s.model_name)
        rsp.fileIO.ensure_dir(match_dict['plots'] + '/')

        match_bg = MatchUtils.make_match_bg_cmd(s.ast_mag1, s.ast_mag2,
                                                outfile=match_dict['bg'])
        print match_bg
        phot = MatchUtils.make_match_bg_cmd(g.mag1, g.mag2,
                                            outfile=match_dict['phot'])

        match_kwargs = {'dmod': s.data.get_col('m - M0')[0],
                        'Av': s.data.get_col('Av')[0],
                        'filter1': s.filter1.replace('F', 'IR'),
                        'filter2': s.filter2.replace('F', 'IR'),
                        'color': s.ast_color,
                        'mag': s.ast_mag2,
                        'pmfile': match_dict['pars']}

        pm_file = MatchUtils.make_calcsfh_param_file(match_bg, **match_kwargs)
        # run match

        match_out = MatchUtils.call_match(pm_file, match_dict['phot'],
                                          s.file_ast, match_dict['output'],
                                          match_dict['msgs'])
        chi, fit = MatchUtils.get_fit(match_out)
        # make plot
        alabel = '$\mathrm{%s}' % s.model_name.replace('_', '\ ')
        metallicity = g.z
        if not np.isfinite(g.z):
            metallicity = '...'
            alabel += '\ Z = %s$' % metallicity
        else:
            alabel += '\ Z = %.4f$' % metallicity
        cmdgrid = match_out + '.cmd'
        labs = ['$%s$' % s.ID, alabel, '$\mathrm{Diff}$',
                '$\chi^2=%g\ \mathrm{Sig}$' % chi]
        grid = MatchUtils.pgcmd(cmdgrid, labels=labs)
        grid[2].set_ylabel(r'$%s$' % g.filter2, fontsize=20)
        grid[2].set_xlabel(r'$%s - %s$' % (g.filter1, g.filter2), fontsize=20)
        figname = os.path.join(match_dict['plots'], os.path.split(cmdgrid)[1])
        print figname
        plt.savefig('%s.png' % figname)
        plt.close()
        s.chi = chi
        s.fit = fit
    return Gals, SGals

def short_list():
    IDs = ['DDO82', 'NGC2403-HALO-6', 'NGC2976-DEEP', 'NGC4163',
           'NGC7793-HALO-6', 'UGC8508', 'UGCA292']
    return IDs

if __name__ == "__main__":

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    logger.info('start of run')
    
    parser = OptionParser()

    usage="%prog model [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-s", action="store_true", default=False,
                      help="run on short set of galaxies")

    parser.add_option("-m", action="store_true", default=False,
                      help="make plots")

    parser.add_option("-p", action="store_true", default=False,
                      help="publish plots")
    
    parser.add_option("-o", action="store_true", default=False,
                      help="out directory set to SNAP/models")

    (options, args) = parser.parse_args()
    
    if not options.s:
        IDs = all_IDs()
    else:
        IDs = short_list()
    
    models = [args[0]]
    
    kwargs = {}
    if options.m:
        kwargs['make_plots'] = True
    if options.p:
        kwargs['publish_plots'] = True
    if options.o:
        kwargs['outdir'] = '/Users/phil/research/TP-AGBcalib/SNAP/'
    
    main(IDs=IDs, models=models, **kwargs)

    '''
    #IDs = [IDs[0]]
    #model = "cmd_input_CAF09_S_SCS.dat"

    models = ["cmd_input_CAF09_s_aug12.dat",
              "cmd_input_CAF09_S_SCS.dat",
              "cmd_input_CAF09_S_SCSFG.dat",
              "cmd_input_CAF09_S_SCSFG_ETA2.dat"]
    #for model in models:
    #    match_tests(IDs, model)
    #    compare_sims(IDs, model)
    chi2_plot(IDs, models)
    sfr_dir = os.path.join(data_src, 'sfh - 0.3dex')
    out_dir = os.path.join(model_src, '0.3dex')
    #plt_dir ...
    kwargs = {'make_plots': True,
              'publish_plots': True,
              'sfr_dir': sfr_dir,
              'out_dir': out_dir}
    main(IDs=IDs, **kwargs)
    '''

    '''
    for model in models:
        sgal, gal = compare_sims(IDs, model)
        fig_loc = os.path.join(plt_dir, sgal.mix, sgal.model_name)
        diag_loc = os.path.join(fig_loc, 'diag')
        html_file = os.path.join(diag_loc, sgal.model_name + '_diag.html')
        gst.side_by_side(diag_loc, sgal.model_name, html_file)
        html_file = os.path.join(fig_loc, sgal.model_name + '.html')
        gst.one_col(fig_loc, sgal.model_name, html_file)
    '''
