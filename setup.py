import setuptools
import os
TOKEN_VALUE = os.getenv('EXPORTED_VAR_WITH_TOKEN')

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
    install_requires=['numpy', 'pandas', 'scipy', 'bottleneck', "Click",
                      f"epiread-tools @ git+https://{TOKEN_VALUE}@github.com/methylgrammarlab/epiread-tools.git"
                      ],
    include_package_data=True,
    entry_points={
    "console_scripts":[
    "deconvolution = deconvolution_models.main:main", #TODO: change to main
    ]
    },
)