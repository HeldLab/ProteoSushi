import codecs
import setuptools

#with open("README.md", "r") as fh:
with codecs.open("README.md", mode='r', encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="proteosushi",
    version="1.2.0",
    author="Rob Seymour",
    author_email="rseymour@wustl.edu",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HeldLab/ProteoSushi",
    download_url="https://github.com/HeldLab/ProteoSushi/archive/v1.1.0.tar.gz",
    include_package_data=True,
    packages=setuptools.find_packages(),
    package_dir={
        "lib": "proteosushi",
        #"": "lib"
    },
    package_data={
        "lib": ['spec_list_fixed.tsv']
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: X11 Applications :: Qt",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "PyQt5>=5.15",
        "requests>=2.24",
        "pandas>=1.0"
    ]
)