from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = 'Gene_POCKET',
    version = '0.0.4',
    keywords='candidate genes prioritization',
    description = 'a library for prioritizing the candidate genes by incorporating information of Knowledge-based gene sets, Effects of variants, GWAS and TWAS',
    license = 'MIT License',
    url = 'https://github.com/zhaouu/POCKET',
    author = 'Hu Zhao',
    author_email = 'zhaohu@webmail.hzau.edu.cn',
    packages = find_packages(),
    include_package_data = True,
    platforms = 'any',
    long_description= long_description,
    long_description_content_type="text/markdown",
    #python_requires='!=3.0.*, !=3.1.*, !=3.2.*,!=3.3.*, !=3.4.*, !=3.5.*, >=3.6.*, <4',
    install_requires=['numpy', 'scipy', 'scikit-learn', 'limix', 'pandas','joblib']
)
