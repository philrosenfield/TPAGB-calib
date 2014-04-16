from .. import tables
import numpy as np
from .. import TPAGBparams

class FileIO(object):
    def __init__(self):
        pass

    def check_target(self, target):
        if target is None:
            assert hasattr(self, 'target'), \
                'need to pass target or have attribute'
        else:
            self.target = target
        return

    def load_data_for_normalization(self, target=None, ags=None):
        self.check_target(target)
        '''load the numbers of data rgb and agb stars from self.ags'''
        if ags is None:
            ags = self.ags

        if len(ags.offsets) == 0:
            self.comp_data = tables.read_completeness_table()

        column_names = ags.data.dtype.names
        if '404' in self.target:
            target = self.target.replace('-deep', '')
        else:
            target = self.target
        row, = np.nonzero(ags.data['target'] == target)
        try:
            [self.__setattr__('%s' % c, ags.data[row]['%s' % c][0])
             for c in column_names if c != 'target']
        except IndexError:
            '%s not found in data table' % target
        if len(ags.offsets) == 0:
            self.ir_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', target.upper(),
                                                      'ir_filter2')
            self.opt_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', target.upper(),
                                                      'opt_filter2')
        else:
            self.ir_offset = self.ir_trgb + ags.offsets[1]
            self.opt_offset = self.opt_trgb + ags.offsets[0]

        self.opt_bins = np.arange(self.opt_min, self.opt_max, 0.1)
        self.ir_bins = np.arange(self.ir_min, self.ir_max, 0.1)

        return self.nopt_rgb, self.nir_rgb

    def load_lf_file(self, lf_file):
        with open(lf_file, 'r') as lff:
            lines = [l.strip() for l in lff.readlines()
                     if not l.startswith('#')]

        hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
        binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]
        return hists, binss

    def load_galaxies(self, hist_it_up=True, target=None, ags=None,
                      color_cut=False, completeness_correction=False):
        self.check_target(target)
        if not hasattr(self, 'opt_bins'):
            self.load_data_for_normalization(ags=ags)

        ir_gal = galaxy_tests.load_galaxy(self.target, band='ir')
        fits_src = TPAGBparams.snap_src + '/data/angst_no_trim'
        opt_gal = galaxy_tests.load_galaxy(self.target, band='opt',
                                           fits_src=TPAGBparams.fits_src)
        # make galaxy histograms
        if color_cut is True:
            filter1 = get_filter1(target)
            opt_color_inds = np.nonzero(opt_gal.color > get_color_cut(filter1))
            ir_color_inds = np.nonzero(ir_gal.color > get_color_cut('F110W'))
            opt_gal.color_cut = opt_color_inds
            ir_gal.color_cut = ir_color_inds
        if hist_it_up is True:
            opt_gal.hist, opt_gal.bins = galaxy_tests.hist_it_up(opt_gal.mag2,
                                                        threash=5)
            ir_gal.hist, ir_gal.bins = galaxy_tests.hist_it_up(ir_gal.mag2,
                                                               threash=5)
        else:
            #nbins = (np.max(opt_gal.mag2) - np.min(opt_gal.mag2)) / 0.1
            if color_cut is False:
                opt_gal.hist, opt_gal.bins = np.histogram(opt_gal.mag2,
                                                          bins=self.opt_bins)
                ir_gal.hist, ir_gal.bins = np.histogram(ir_gal.mag2,
                                                        bins=self.ir_bins)
            else:
                opt_gal.hist, opt_gal.bins = \
                    np.histogram(opt_gal.mag2[opt_color_inds],
                                 bins=self.opt_bins)
                ir_gal.hist, ir_gal.bins = \
                    np.histogram(ir_gal.mag2[ir_color_inds], bins=self.ir_bins)

        if completeness_correction is True:
            ast_dict = tables.read_completeness_corrections()
            opt_correction = ast_dict[target]['opt_correction'][:-1]
            ir_correction = ast_dict[target]['ir_correction'][:-1]
            opt_gal.hist *= 1. / opt_correction
            ir_gal.hist *= 1. / ir_correction

        return opt_gal, ir_gal

    def load_trilegal_data(self):
        '''load trilegal F814W and F160W mags'''

        if hasattr(self, 'target'):
            filter1 = get_filter1(self.target)
        else:
            print 'help, I need filter1!!'
        opt_mag = self.sgal.data.get_col('F814W')
        ir_mag = self.sgal.data.get_col('F160W')
        opt_mag1 = self.sgal.mag1
        ir_mag1 = self.sgal.data.get_col('F110W')

        opt_color_cut, = \
            np.nonzero((opt_mag1 - opt_mag) > get_color_cut(filter1))
        ir_color_cut, = \
            np.nonzero((ir_mag1 - ir_mag) > get_color_cut('F110W'))
        self.opt_color_cut = opt_color_cut
        self.ir_color_cut = ir_color_cut
        self.shift_mags(opt_inds=opt_color_cut, ir_inds=ir_color_cut)
        return

    def read_trilegal_catalog(self, trilegal_output, filter1='F606W'):
        '''read the trilegal cat mag1 and mag2 are optical.'''
        self.sgal = rsp.Galaxies.simgalaxy(trilegal_output, filter1=filter1,
                                           filter2='F814W')
        return
