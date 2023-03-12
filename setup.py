from pathlib import Path
from setuptools import find_packages, setup
import versioneer

PACKAGE_DIR = "."
install_requires = Path("requirements.txt").read_text().split("\n")

setup(
    package_dir = {"": PACKAGE_DIR},
    install_requires = install_requires,
    packages = find_packages(where=PACKAGE_DIR),
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)