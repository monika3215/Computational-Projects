import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1

adql_query = '''
select top 2000 parallax, phot_g_mean_mag, bp_rp, parallax_over_error ,ruwe
from gaiadr3.gaia_source
where parallax>10
  and parallax_over_error>10
  and ruwe <1.4
  and bp_rp is not null'''

job = Gaia.launch_job(adql_query)
r = job.get_results().to_pandas()

G= r['phot_g_mean_mag']
parallax = r['parallax']
bp_rp=r['bp_rp']
# M_G is the absolute magnitude in the G band
distance_pc = 1000 / parallax
M_G = G + 5*np.log10 (parallax) -10
print(r.head(6))
print(len(r))
color = bp_rp 

plt.figure(figsize=(8, 6))
plt.scatter(color,M_G,c='red',s=3,alpha=0.5)
plt.gca().invert_yaxis()
plt.xlabel("Color Index($bp-rp$)")
plt.ylabel("Absolute Magnitude ($M_G$)")
plt.title("HR Diagram of Stars (Gaia DR3)")
plt.show()