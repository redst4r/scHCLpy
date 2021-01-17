from setuptools import setup, find_packages

setup(name='scHCLpy',
      version=0.1,
      description='A simplistic python version of the scHCL R package',
      url='http://github.com/redst4r/scHCLpy/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='scanpy, scrnaseq',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          ],
      zip_safe=False)
