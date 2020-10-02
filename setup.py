from setuptools import setup, find_packages

setup(
    name = 'pockt',
    version = '0.0.1',
    keywords='candidate genes prioritization',
    description = 'a library for prioritizing the candidate genes by incorporating information of Knowledge-based gene sets, Effects of variants, GWAS and TWAS',
    license = 'MIT License',
    url = 'https://github.com/zhaouu/POCKT',
    author = 'Hu Zhao',
    author_email = 'zhaohu@webmail.hzau.edu.cn',
    packages = find_packages(),
    include_package_data = True,
    platforms = 'any',
    #python_requires='!=3.0.*, !=3.1.*, !=3.2.*,!=3.3.*, !=3.4.*, !=3.5.*, >=3.6.*, <4',
    install_requires=['numpy', 'scipy', 'scikit-learn', 'limix', 'pandas']
)
