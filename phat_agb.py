import os
import ResolvedStellarPops as rsp
import sfh_tests
# convert MATCH SFH to Trilegal AMR
sfh_loc = '/home/phil/research/TP-AGBcalib/PHAT/sfh'
sfh_files = rsp.fileIO.get_files(sfh_loc, '*sfh')

sfhs = [sfh_tests.StarFormationHistories(sfh_file, 'match-old') for sfh_file in sfh_files]
tri_sfhs = [sfh.make_trilegal_sfh() for sfh in sfhs]
trilegal_outputs = []

# values taken from .sfh file!
dmod = 24.47
dist = 10 ** (dmod / 5 + 1.)
Av = 0.4
gal_inp_dict = {'mag_limit_val': 28, 'object_dist': dist, 'object_av': Av,
                'file_imf': 'tab_imf/imf_kroupa_orig.dat'}

cmd_input_files = ['cmd_input_CAF09_S_NOV13.dat', 'cmd_input_CAF09_S_OCT13.dat']
for cmd_input_file in cmd_input_files:
    agb_mod = cmd_input_file.split('_')[-1].replace('.dat', '')
    for tri_sfh in tri_sfhs:
        for aringer in [True, False]:
            # make galaxy input for trilegal
            gal_file = tri_sfh.replace('tri.dat', '%s.gal.inp' % agb_mod.lower())
            if aringer is True:
                gal_file = gal_file.replace('gal.inp', 'aringer.gal.inp')
            # initialze the dict
            gal_init = {'photsys': 'phat_agb', 'filter1': 'F475W', 'object_sfr_file': tri_sfh,
                        'aringer': aringer, 'object_mass': 1e8}
            gal_inp = rsp.fileIO.input_parameters(default_dict=rsp.TrilegalUtils.galaxy_input_dict(**gal_init))
            gal_inp.add_params(gal_inp_dict)
            # write the file
            gal_inp.write_params(gal_file, rsp.TrilegalUtils.galaxy_input_fmt())
            # run trilegal
            trilegal_output = gal_file.replace('gal.inp', 'tri.out')
            if not os.path.isfile(trilegal_output):
                if not os.path.isfile(trilegal_output + '.gz'):
                    rsp.TrilegalUtils.run_trilegal(cmd_input_file, gal_file,
                                               trilegal_output, rmfiles=False,
                                               dry_run=False)
            trilegal_outputs.append(trilegal_output)
print trilegal_outputs