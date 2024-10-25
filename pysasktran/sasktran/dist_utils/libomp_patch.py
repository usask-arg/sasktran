from pathlib import Path
import sys
import pkgutil


def patch_libomp():
    print("Trying to patch openmp")
    # Only on mac
    if sys.platform != 'darwin':
        print("Not on mac, skipping")
        return

    # Check if we have bundled .dylibs
    package = pkgutil.get_loader('sasktran')

    bundled_dylib = Path(package.get_filename()).parent.joinpath('.dylibs/libomp.dylib')

    # First we need to determine if we are in an anaconda environment, and that if libomp exists
    conda_dylib = Path(sys.base_prefix).joinpath('lib/libomp.dylib')

    print(f"conda_dylib: {conda_dylib}")
    print(f"bundled_dylib: {bundled_dylib}")

    if conda_dylib.exists() and bundled_dylib.exists() and not (bundled_dylib.is_symlink()):
        print("Patching libomp.dylib")
        bundled_dylib.unlink()
        bundled_dylib.symlink_to(conda_dylib)
