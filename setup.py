from setuptools import setup, find_packages

setup(
    name='Polyploid_population_simulation',
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
    py_modules=["simulation_pop"],
    author='Z C',
    author_email='cjun01@gmail.com',
    description='A library to simulate and visualize population data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/cjun01/Polyploid_population_simulation', # link to your github or wherever you host the code.
)
