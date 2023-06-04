'''
MIT License

Copyright (c) 2022 irene unterman and ben berman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import setuptools
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='deconvolution_models',
    version='0.0.2',
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
    install_requires=['numpy', 'pandas', 'scipy', 'bottleneck', "Click", "pytest",
                      f"epiread-tools @ git+https://github.com/methylgrammarlab/epiread-tools.git",
                      ],
    include_package_data=True,
    entry_points={
    "console_scripts":[
    "deconvolution = deconvolution_models.main:main",
    ]
    },
)
# python3 setup.py sdist bdist_wheel
