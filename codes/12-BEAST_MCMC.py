import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Updated to the correct log file name
log_file = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_beast_full_v27.log"

print("--- Loading BEAST MCMC Log ---")
# 1. Load the data, ignoring BEAST's header comments
df = pd.read_csv(log_file, sep='\t', comment='#')

# 2. Remove 10% Burn-in
burnin = int(len(df) * 0.1)
df_post = df.iloc[burnin:].copy()

print(f"Total states recorded: {len(df)}")
print(f"States analyzed after 10% burn-in: {len(df_post)}")
print("-" * 40)

# 3. Fast ESS Calculator
def calculate_ess(x):
    n = len(x)
    mu = np.mean(x)
    var = np.var(x, ddof=1)
    if var == 0: return np.nan

    rho_sum = 0.0
    for lag in range(1, n):
        # Calculate autocorrelation at current lag
        cov = np.sum((x[:-lag] - mu) * (x[lag:] - mu)) / (n - 1)
        rho = cov / var
        if rho < 0.05: # Stop summing when correlation drops to noise
            break
        rho_sum += rho

    ess = n / (1.0 + 2.0 * rho_sum)
    return ess

# 4. Check key parameters
params_to_check = ['posterior', 'TreeHeight.t:alignment', 'clockRate.c:alignment']

for p in params_to_check:
    if p in df_post.columns:
        mean_val = df_post[p].mean()
        ess_val = calculate_ess(df_post[p].values)

        status = "✅ GOOD (>200)" if ess_val > 200 else "⚠️ LOW (Run longer or ignore if secondary)"

        print(f"Parameter: {p}")
        print(f"  Mean Value: {mean_val:.4e}")
        print(f"  ESS Score:  {ess_val:.1f}  {status}")
        print("-" * 40)

# 5. Plot the Traces (The "Tracer" View)
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

if 'TreeHeight.t:alignment' in df_post.columns:
    axes[0].plot(df_post['Sample'], df_post['TreeHeight.t:alignment'], color='#1f77b4', alpha=0.8, linewidth=0.5)
    axes[0].set_title('MCMC Trace: TreeHeight (Age of Root Ancestor)')
    axes[0].set_xlabel('State')
    axes[0].set_ylabel('Years')
    axes[0].grid(True, linestyle='--', alpha=0.5)

if 'clockRate.c:alignment' in df_post.columns:
    axes[1].plot(df_post['Sample'], df_post['clockRate.c:alignment'], color='#ff7f0e', alpha=0.8, linewidth=0.5)
    axes[1].set_title('MCMC Trace: Clock Rate')
    axes[1].set_xlabel('State')
    axes[1].set_ylabel('Mutations / site / year')
    axes[1].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()
