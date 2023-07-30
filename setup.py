from setuptools import setup, find_packages

setup(
    name='population_simulation',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'seaborn',
        'tqdm',
        'numba',
        'scikit-learn',
    ],
    author='Z C',
    author_email='cjun01@gmail.com',
    description='A library to simulate and visualize population data.',
    long_description=open('README.md').read(),
    license='MIT'

)
