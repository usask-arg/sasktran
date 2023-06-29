_old_print_mode = None


def pytest_configure():
    global _old_print_mode
    import numpy as np
    from pkg_resources import parse_version

    if parse_version(np.__version__) >= parse_version('1.14.0'):
        # Numpy default print repr changed, breaking some doctests
        OLD_PRINT_MODE = np.get_printoptions()['legacy']
        np.set_printoptions(legacy='1.13')


def pytest_unconfigure():
    global _old_print_mode
    import numpy as np
    from pkg_resources import parse_version

    if parse_version(np.__version__) >= parse_version('1.14.0'):
        # Numpy default print repr changed, breaking some doctests
        np.set_printoptions(legacy=_old_print_mode)