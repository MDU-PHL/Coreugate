from setuptools import setup, find_packages

setup (name = 'COREugate',
       version = '1.0',
       include_package_data = True,
       packages=find_packages(),
       description = 'calculate pair-wise allelic distances from cgMLST',
       author = 'Kristy Horan',
       url = 'https://github.com/kristyhoran/Coreugate',
       install_requires = ['jinja2','biopython>=1.70','pandas>=0.23.0', 'pathlib', 'cleo','numpy', 'snakemake>=5.3.0'],
       python_requires='>=3.6',
       entry_points={
        "console_scripts": [
            "coreugate = run_coreugate:main"
        ],},
        setup_requires=['pytest-runner'],
        tests_require = ['pytest'],
        test_suite = 'test'
)