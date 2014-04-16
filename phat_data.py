
def prepare_vsfh_run(targets, cmd_input_files, nsfhs, mk_tri_sfh_kw=None,
                     vary_sfh_kw=None, make_many_kw=None, vsfh_kw=None,
                     table_file='default', dry_run=False, default_kw=None):
    '''
    Run a number of SFH variations on a galaxy.
    If passed default to the args, will attempt to find the file based on the
    galaxy name.

    ARGS:
    galaxy_name: target name, ex: ddo71 (case doesn't matter)
    cmd_input_file: filename, ex: 'cmd_input_CAF09_S_OCT13.dat'
    match_sfh_file: 'default', the sfh file from match.
    match_fileorigin: 'match-grid', which type of sfh file
                        'match' or 'match-grid'
    galaxy_input_file: 'default', base input file for trilegal, will be copied.

    mk_tri_sfh_kw: A dict to be passed to VarySFH.make_trilegal_sfh
        default: random_sfh = True, random_z = False

    make_many_kw: A dict to be passed to VarySFH.prepare_trilegal_sfr
        default: nsfhs = 50, mk_tri_sfh_kw dict.

    vary_sfh_kw: A dict to be passed to VarySFH.vary_the_sfh
        default: diag_plots = True, make_many_kw dict

    RETURNS:
    VarySFHs class
    '''
    # load and possibly overwrite make_many_kw from defaults
    # start the stream logger. The file logger is done target by target.
    global logger
    logger = logging.getLogger()
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S",
                                                   time.localtime()))
    # create formatter and add it to the handlers
    formatter = \
        '%(asctime)-15s %(levelname)s %(funcName)s %(lineno)d %(message)s'
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S",
                                                   time.localtime()))

    make_many_kw = make_many_kw or {}
    vary_sfh_kw = vary_sfh_kw or {}
    mk_tri_sfh_kw = mk_tri_sfh_kw or {}
    default_kw = default_kw or {}

    mk_tri_sfh_kw = dict({'dry_run': dry_run}.items() + mk_tri_sfh_kw.items())

    make_many_kw = dict({'mk_tri_sfh_kw': mk_tri_sfh_kw}.items()
                         + make_many_kw.items())

    vary_sfh_kw = dict({'make_many_kw': make_many_kw}.items()
                        + vary_sfh_kw.items())

    vsfh_kw = vsfh_kw or {}
    vsfh_kw = dict({'file_origin': 'match-hmc', 'table_file': table_file,
                    'outfile_loc': 'default', 'nsfhs': nsfhs}.items() \
                   + vsfh_kw.items())

    vsfh_kws = []
    vSFHs = []
    for target, cmd_input in itertools.product(targets, cmd_input_files):
        target = target.lower()
        #print target, cmd_input
        vsfh_kws.append(default_values_for_vsfh(target, cmd_input,
                                                vsfh_kw=vsfh_kw, **default_kw))
        vSFHs.append(VarySFHs(**vsfh_kw))
    return vSFHs, vary_sfh_kw