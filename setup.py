import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='deconvolution_models',
    version='0.0.1',
    author='Irene Unterman',
    author_email='irene.guberman@mail.huji.ac.il',
    description='Models for WGBS deconvolution',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/methylgrammarlab/deconvolution_models',
    project_urls = {
        "Bug Tracker": "https://github.com/methylgrammarlab/deconvolution_models/issues"
    },
    license='MIT',
    packages=['deconvolution_models'],
    install_requires=['numpy', 'pandas', 'scipy', 'bottleneck', "Click"],
    include_package_data=True,
    entry_points={
    "console_scripts":[
    "celfie = deconvolution_models.celfie:main",
    "celfie-plus = deconvolution_models.celfie_plus:main",
    "epistate = deconvolution_models.epistate:main",
    "epistate-plus = deconvolution_models.epistate_plus:main",
    ]
    },
)