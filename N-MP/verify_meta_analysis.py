"""
Meta-analysis verification script
Random-effects model (DerSimonian-Laird) with Freeman-Tukey double arcsine transformation
"""
import numpy as np
from scipy import stats

studies = [
    {'id':1,  'author':'Leonard et al.',       'year':2024, 'events':18, 'n':20,  'tissue':'Blood/Serum/Plasma', 'pmid':'38761430', 'country':'UK'},
    {'id':2,  'author':'Ke et al.',             'year':2023, 'events':59, 'n':69,  'tissue':'Stool/Gut',          'pmid':'37837933', 'country':'China'},
    {'id':3,  'author':'Arslan et al.',         'year':2025, 'events':10, 'n':12,  'tissue':'Cerumen/Ear',        'pmid':'40933963', 'country':'Turkey'},
    {'id':4,  'author':'Santini et al.',        'year':2024, 'events':5,  'n':6,   'tissue':'Stool/Gut',          'pmid':'38883588', 'country':'Italy'},
    {'id':5,  'author':'Zakynthinos et al.',    'year':2025, 'events':5,  'n':8,   'tissue':'Lung',               'pmid':'41450414', 'country':'Greece'},
    {'id':6,  'author':'Demirelli et al.',      'year':2024, 'events':6,  'n':12,  'tissue':'Semen/Testis',       'pmid':'38745203', 'country':'Turkey'},
    {'id':7,  'author':'Saraluck et al.',       'year':2024, 'events':23, 'n':59,  'tissue':'Breast milk',        'pmid':'39064070', 'country':'Thailand'},
    {'id':8,  'author':'Ozgen Alpaydin et al.', 'year':2024, 'events':10, 'n':18,  'tissue':'Lung',               'pmid':'39368622', 'country':'Turkey'},
    {'id':9,  'author':'Halfar et al.',         'year':2023, 'events':9,  'n':10,  'tissue':'Placenta',           'pmid':'37769922', 'country':'Czech Republic'},
    {'id':10, 'author':'Qu et al.',             'year':2025, 'events':111,'n':200, 'tissue':'Semen/Testis',       'pmid':'41068921', 'country':'China'},
    {'id':11, 'author':'Padarya et al.',        'year':2025, 'events':42, 'n':60,  'tissue':'Lung',               'pmid':'41523016', 'country':'India'},
]

k = len(studies)
N_total = sum(s['n'] for s in studies)

def ft_transform(x, n):
    """Freeman-Tukey double arcsine transformation"""
    return np.arcsin(np.sqrt(x / (n + 1))) + np.arcsin(np.sqrt((x + 1) / (n + 1)))

def ft_variance(n):
    """Variance of FT-transformed proportion"""
    return 1.0 / (n + 0.5)

def ft_backtransform(t):
    """Back-transform FT to proportion"""
    return (np.sin(t / 2)) ** 2

# Transform all studies
yi = np.array([ft_transform(s['events'], s['n']) for s in studies])
vi = np.array([ft_variance(s['n']) for s in studies])
ni = np.array([s['n'] for s in studies])

# Fixed-effects
wi_fe = 1.0 / vi
theta_fe = np.sum(wi_fe * yi) / np.sum(wi_fe)

# Cochran's Q
Q = np.sum(wi_fe * (yi - theta_fe) ** 2)
df = k - 1
Q_pval = 1.0 - stats.chi2.cdf(Q, df)
I2 = max(0, (Q - df) / Q) * 100

# DerSimonian-Laird tau^2
C = np.sum(wi_fe) - np.sum(wi_fe ** 2) / np.sum(wi_fe)
tau2 = max(0, (Q - df) / C)

# Random-effects
wi_re = 1.0 / (vi + tau2)
theta_re = np.sum(wi_re * yi) / np.sum(wi_re)
se_re = np.sqrt(1.0 / np.sum(wi_re))

# Back-transform
pooled_p = ft_backtransform(theta_re)
ci_lo = ft_backtransform(theta_re - 1.96 * se_re)
ci_hi = ft_backtransform(theta_re + 1.96 * se_re)

# Prediction interval
t_crit = stats.t.ppf(0.975, df)
pi_lo = ft_backtransform(theta_re - t_crit * np.sqrt(tau2 + se_re ** 2))
pi_hi = ft_backtransform(theta_re + t_crit * np.sqrt(tau2 + se_re ** 2))

print("=" * 70)
print("META-ANALYSIS VERIFICATION RESULTS")
print("=" * 70)
print(f"Studies (k): {k}")
print(f"Total participants (N): {N_total}")
print()
print("OVERALL POOLED ESTIMATE (Random-Effects, DerSimonian-Laird):")
print(f"  Prevalence: {pooled_p * 100:.1f}%")
print(f"  95% CI: [{ci_lo * 100:.1f}%, {ci_hi * 100:.1f}%]")
print(f"  Prediction Interval: [{pi_lo * 100:.1f}%, {pi_hi * 100:.1f}%]")
print()
print("HETEROGENEITY:")
print(f"  Q = {Q:.2f}, df = {df}, p = {Q_pval:.4f}")
print(f"  I-squared = {I2:.1f}%")
print(f"  tau-squared = {tau2:.4f}")
print()

# Study-level details
print("INDIVIDUAL STUDY DETAILS:")
print(f"{'Author':<30s} {'Year':>4s} {'Events/N':>8s} {'Prev%':>6s} {'Weight%':>7s}")
print("-" * 60)
total_w = np.sum(wi_re)
for i, s in enumerate(studies):
    p_raw = s['events'] / s['n'] * 100
    w_pct = wi_re[i] / total_w * 100
    print(f"{s['author']:<30s} {s['year']:>4d} {s['events']:>3d}/{s['n']:<3d}  {p_raw:>5.1f}% {w_pct:>6.1f}%")

# Egger's test
print()
print("PUBLICATION BIAS (Egger's Test):")
sei = np.sqrt(vi)
# Standard Egger: regress t_i = y_i/se_i on 1/se_i
ti_egger = yi / sei
precision = 1.0 / sei
slope, intercept, r, p_egger, se_slope = stats.linregress(precision, ti_egger)
print(f"  Intercept: {intercept:.3f}")
print(f"  p-value: {p_egger:.4f}")
print(f"  Conclusion: {'Publication bias detected' if p_egger < 0.05 else 'No significant publication bias'}")

# Leave-one-out
print()
print("LEAVE-ONE-OUT SENSITIVITY ANALYSIS:")
print(f"{'Excluded Study':<30s} {'Pooled%':>7s} {'95% CI':>18s} {'I2':>5s}")
print("-" * 65)
for excl in range(k):
    mask = np.ones(k, dtype=bool)
    mask[excl] = False
    yi_loo = yi[mask]
    vi_loo = vi[mask]
    k_loo = k - 1
    wi_loo = 1.0 / vi_loo
    theta_fe_loo = np.sum(wi_loo * yi_loo) / np.sum(wi_loo)
    Q_loo = np.sum(wi_loo * (yi_loo - theta_fe_loo) ** 2)
    I2_loo = max(0, (Q_loo - (k_loo - 1)) / Q_loo) * 100 if Q_loo > 0 else 0
    C_loo = np.sum(wi_loo) - np.sum(wi_loo ** 2) / np.sum(wi_loo)
    tau2_loo = max(0, (Q_loo - (k_loo - 1)) / C_loo)
    wi_re_loo = 1.0 / (vi_loo + tau2_loo)
    theta_re_loo = np.sum(wi_re_loo * yi_loo) / np.sum(wi_re_loo)
    se_re_loo = np.sqrt(1.0 / np.sum(wi_re_loo))
    p_loo = ft_backtransform(theta_re_loo) * 100
    ci_l = ft_backtransform(theta_re_loo - 1.96 * se_re_loo) * 100
    ci_h = ft_backtransform(theta_re_loo + 1.96 * se_re_loo) * 100
    print(f"{studies[excl]['author']:<30s} {p_loo:>6.1f}% [{ci_l:>5.1f}%, {ci_h:>5.1f}%] {I2_loo:>4.0f}%")

# Subgroup by tissue
print()
print("SUBGROUP ANALYSIS BY TISSUE:")
tissues = {}
for i, s in enumerate(studies):
    t = s['tissue']
    if t not in tissues:
        tissues[t] = []
    tissues[t].append(i)

print(f"{'Tissue':<25s} {'Pooled%':>7s} {'95% CI':>18s} {'k':>3s} {'I2':>5s}")
print("-" * 65)
for tissue, indices in sorted(tissues.items()):
    k_sub = len(indices)
    if k_sub == 1:
        idx = indices[0]
        x, n = studies[idx]['events'], studies[idx]['n']
        p_hat = x / n
        z = 1.96
        denom = 1 + z ** 2 / n
        center = (p_hat + z ** 2 / (2 * n)) / denom
        margin = z * np.sqrt(p_hat * (1 - p_hat) / n + z ** 2 / (4 * n ** 2)) / denom
        ci_l = max(0, center - margin) * 100
        ci_h = min(1, center + margin) * 100
        print(f"{tissue:<25s} {p_hat * 100:>6.1f}% [{ci_l:>5.1f}%, {ci_h:>5.1f}%] {k_sub:>3d}   0%")
    else:
        yi_sub = yi[np.array(indices)]
        vi_sub = vi[np.array(indices)]
        wi_sub = 1.0 / vi_sub
        theta_fe_sub = np.sum(wi_sub * yi_sub) / np.sum(wi_sub)
        Q_sub = np.sum(wi_sub * (yi_sub - theta_fe_sub) ** 2)
        df_sub = k_sub - 1
        I2_sub = max(0, (Q_sub - df_sub) / Q_sub) * 100 if Q_sub > 0 else 0
        C_sub = np.sum(wi_sub) - np.sum(wi_sub ** 2) / np.sum(wi_sub)
        tau2_sub = max(0, (Q_sub - df_sub) / C_sub)
        wi_re_sub = 1.0 / (vi_sub + tau2_sub)
        theta_re_sub = np.sum(wi_re_sub * yi_sub) / np.sum(wi_re_sub)
        se_re_sub = np.sqrt(1.0 / np.sum(wi_re_sub))
        p_sub = ft_backtransform(theta_re_sub) * 100
        ci_l = ft_backtransform(theta_re_sub - 1.96 * se_re_sub) * 100
        ci_h = ft_backtransform(theta_re_sub + 1.96 * se_re_sub) * 100
        print(f"{tissue:<25s} {p_sub:>6.1f}% [{ci_l:>5.1f}%, {ci_h:>5.1f}%] {k_sub:>3d} {I2_sub:>4.0f}%")

print()
print("=" * 70)
print("CROSS-VALIDATION WITH ORIGINAL REPORT:")
print(f"  Original pooled: 68.5% [56.3-78.6%], I2=77%")
print(f"  Verified pooled: {pooled_p * 100:.1f}% [{ci_lo * 100:.1f}-{ci_hi * 100:.1f}%], I2={I2:.0f}%")
print(f"  MATCH: {'YES' if abs(pooled_p * 100 - 68.5) < 2.0 else 'DISCREPANCY - CHECK'}")
print("=" * 70)
