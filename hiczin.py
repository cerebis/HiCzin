import pandas
import numpy as np
import patsy
import argparse
import statsmodels.discrete.count_model as cm
import seaborn as sb

from scipy.stats import nbinom

#
# TODO -- this should really be part of bin3C or use proxigenomics_toolkit
#         and derive all this information itself. The logic for producing the
#         contacts is important. An existing issue is that bin3C is Python 2.7
#         while I strongly suspect the latest codebase of important modules
#         will only be deployed for Python >3
#

def fit_data(df, n_sample=None):
    if n_sample is None:
        return df
    else:
        return df.sample(n=n_sample)


def mul_ab(df, index1, index2, col):
    return df.loc[index1, col].values * df.loc[index2, col].values


def scaler(a, _mu=None, _sigma=None):
    if not _mu:
        _mu = a.mean()
    if not _sigma:
        _sigma = a.std()
    return (a - float(_mu)) / float(_sigma), _mu, _sigma


def make_table(_contacts, _sites, _lengths, _coverage):
    return pandas.DataFrame({'contacts': _contacts,
                             'sites': _sites,
                             'length': _lengths,
                             'coverage': _coverage})


def convert_params(mu, alpha):
    """
    Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports

    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
    """
    r = 1. / alpha
    var = mu + 1. / r * mu ** 2
    p = (var - mu) / var
    return r, 1. - p


def nbinom_cdf(c, mu, alpha):
    """
    Re-parameterised scipy negative binomial

    :param c: observed count
    :param mu: expected count
    :param alpha: dispersion (alpha = 1/r)
    :return: cumulative probability
    """
    return nbinom.cdf(c, *convert_params(mu, alpha))


parser = argparse.ArgumentParser()
parser.add_argument('--threshold', '-t', type=float, default=0.1,
                    help='Threshold cut-off for significant contacts')
parser.add_argument('--n-sample', '-n', type=int,
                    help='Limit number of observations used in model fitting')
parser.add_argument('coverage', help='Coverage file [name,depth]')
parser.add_argument('contig_info', help="Contig information file [name,length,sites]")
parser.add_argument('intra_contacts', help='Intra-genome contacts')
parser.add_argument('raw_contacts', help='All raw contacts')
parser.add_argument('valid_out', help='Output table of valid contacts')
parser.add_argument('spur_out', help='Output table of spurious contacts')
args = parser.parse_args()

print('Reading input data')

coverage = pandas.read_csv(args.coverage, header=None, names=['contig_name', 'coverage'])
contig_info = pandas.read_csv(args.contig_info, header=None, names=['contig_name', 'length', 'sites'])
sample_data = pandas.read_csv(args.intra_contacts, header=None, names=['index1', 'index2', 'contacts'])
all_data = pandas.read_csv(args.raw_contacts, header=None, names=['index1', 'index2', 'contacts'])

# combine coverage with contig info.
# TODO this should be a single table, if we were to keep it separate
contig_info['sites'] = contig_info.sites.replace(0, 1)
coverage['coverage'] = coverage.coverage.replace(0, 1)
contig_info = contig_info.set_index('contig_name').join(coverage.set_index('contig_name'))

print('Preparing variables')

# data preparation involves transformations to approximate normal
sample_sites = np.log(mul_ab(contig_info, sample_data.index1, sample_data.index2, 'sites'))
sample_length = np.log(mul_ab(contig_info, sample_data.index1, sample_data.index2, 'length'))
sample_coverage = np.sqrt(np.log(mul_ab(contig_info, sample_data.index1, sample_data.index2, 'coverage')))
sample_contacts = sample_data.contacts.values.astype(np.int32)

all_sites = np.log(mul_ab(contig_info, all_data.index1, all_data.index2, 'sites'))
all_length = np.log(mul_ab(contig_info, all_data.index1, all_data.index2, 'length'))
all_coverage = np.sqrt(np.log(mul_ab(contig_info, all_data.index1, all_data.index2, 'coverage')))
all_contacts = all_data.contacts.values.astype(np.int32)

# simple data standardisation
# apply the same transformation to the raw data columns
sample_sites, _mu, _sigma, = scaler(sample_sites)
all_site = scaler(all_sites, _mu, _sigma)[0]

sample_length, _mu, _sigma = scaler(sample_length)
all_length = scaler(all_length, _mu, _sigma)[0]

sample_coverage, _mu, _sigma = scaler(sample_coverage)
all_coverage = scaler(all_coverage, _mu, _sigma)[0]

# lets see if these variables look standard normal
sb.set_theme(style="white")
df = pandas.DataFrame({'length': sample_length, 'sites': sample_sites, 'coverage': sample_coverage}).sample(1000)
g = sb.PairGrid(df, diag_sharey=False)
g.map_upper(sb.scatterplot, alpha=0.5, s=15)
g.map_lower(sb.kdeplot)
g.map_diag(sb.kdeplot, lw=2)
g.savefig('data_distrib.pdf')


sample_exog = make_table(sample_contacts, sample_sites, sample_length, sample_coverage)

print('Fitting ZINB model to data')

# prepare design matrices to supply ZINB model
formula = "contacts ~ sites * length + coverage"
y, X = patsy.dmatrices(formula, fit_data(sample_exog, args.n_sample), return_type='matrix')

# attempt to fit
md = cm.ZeroInflatedNegativeBinomialP(y,  # endogenous
                                      X,  # exogenous
                                      X,  # zi
                                      method='l1', inflation='logit')
mdf = md.fit_regularized()
print(mdf.summary())

print('Calculating results')

# prediction means include zero-inflation adjustment
#  -- this was not part of the HiCzin R code, which ignores ZI part
_, X = patsy.dmatrices(formula, sample_exog, return_type='matrix')
sample_mu = mdf.predict(exog=X, exog_infl=X, which='mean')
sample_norm = sample_contacts / sample_mu
ix_nz = sample_norm > 0
sample_nz_mu = sample_mu[ix_nz]
sample_nz_norm = sample_norm[ix_nz]

all_exog = make_table(all_contacts, all_sites, all_length, all_coverage)
_, X = patsy.dmatrices(formula, all_exog, return_type='matrix')
all_mu = mdf.predict(exog=X, exog_infl=X, which='mean')
all_norm = all_contacts / all_mu

alpha = mdf.params[-1]

sample_pvalue = nbinom_cdf(sample_contacts[ix_nz], sample_nz_mu, alpha)
all_pvalue = nbinom_cdf(all_contacts, all_mu, alpha)

# spurious are marked as low p-value **OR** low normalised contact
# TODO Note: very high observed compared to expected will receive low p-values
#            and be excluded. Are these truly bad?
ix_spur = ((all_pvalue < np.quantile(sample_pvalue, args.threshold)) |
           (all_norm < np.quantile(sample_nz_norm, args.threshold)))

# add new result columns
all_data['expected'] = all_mu
all_data['normed'] = all_norm
all_data['pvalue'] = all_pvalue

# add additional information for each index
all_data[['length1', 'sites1', 'coverage1']] = \
    contig_info.loc[all_data.index1, ['length', 'sites', 'coverage']].reset_index(drop=True)
all_data[['length2', 'sites2', 'coverage2']] = \
    contig_info.loc[all_data.index2, ['length', 'sites', 'coverage']].reset_index(drop=True)

# split raw interactions into accepted and spurious
all_valid = all_data[~ix_spur]
all_spur = all_data[ix_spur]

print('Writing output')

all_valid.to_csv(args.valid_out,
                 index=False, header=True)
all_spur.to_csv(args.spur_out,
                index=False, header=True)
