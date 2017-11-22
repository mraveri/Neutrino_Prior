"""

This file contains an example on how to use the neutrino fitting formula
to produce some prior plots that appeared in the paper.

Developed by: Marco Raveri

"""

import neutrino_fitting_form as nup
import matplotlib.pyplot as plt
import numpy as np

outdir = './results/'
x = np.linspace( 0.0, 0.3, 1000 )

# plot figure 6 panel a:
plt.plot( x, nup.evidence_M_NO/(nup.evidence_M_NO+nup.evidence_M_IO)*nup.fitting_function_M_NO(x), label='Majorana NO' );
plt.plot( x, nup.evidence_D_NO/(nup.evidence_D_NO+nup.evidence_D_IO)*nup.fitting_function_D_NO(x), label='Dirac NO' );
plt.plot( x, nup.evidence_S_NO/(nup.evidence_S_NO+nup.evidence_S_IO)*nup.fitting_function_S_NO(x), label='Seesaw NO' );
plt.xlabel('$\\Sigma m_\\nu \,\,[\\, {\\rm eV}\\,]$');
plt.ylabel('$\\pi(  \\Sigma m_\\nu) \,\,[\\, {\\rm eV}^{-1}\\,]$');
plt.legend();
plt.savefig(outdir+'1_figure6a.pdf')
plt.clf()

# plot figure 6 panel b:
plt.plot( x, nup.evidence_M_IO/(nup.evidence_M_NO+nup.evidence_M_IO)*nup.fitting_function_M_IO(x), label='Majorana IO' );
plt.plot( x, nup.evidence_D_IO/(nup.evidence_D_NO+nup.evidence_D_IO)*nup.fitting_function_D_IO(x), label='Dirac IO' );
plt.plot( x, nup.evidence_S_IO/(nup.evidence_S_NO+nup.evidence_S_IO)*nup.fitting_function_S_IO(x), label='Seesaw IO' );
plt.xlabel('$\\Sigma m_\\nu \,\,[\\, {\\rm eV}\\,]$');
plt.ylabel('$\\pi(  \\Sigma m_\\nu) \,\,[\\, {\\rm eV}^{-1}\\,]$');
plt.legend();

plt.savefig(outdir+'1_figure6b.pdf')
plt.clf()

# plot figure 7:
plt.plot( x, nup.evidence_M_NO/(nup.evidence_M_NO+nup.evidence_M_IO)*nup.fitting_function_M_NO(x)
             +nup.evidence_M_IO/(nup.evidence_M_NO+nup.evidence_M_IO)*nup.fitting_function_M_IO(x), label='Majorana NO+IO' );
plt.plot( x, nup.evidence_D_NO/(nup.evidence_D_NO+nup.evidence_D_IO)*nup.fitting_function_D_NO(x)
             +nup.evidence_D_IO/(nup.evidence_D_NO+nup.evidence_D_IO)*nup.fitting_function_D_IO(x), label='Dirac NO+IO' );
plt.plot( x, nup.evidence_S_NO/(nup.evidence_S_NO+nup.evidence_S_IO)*nup.fitting_function_S_NO(x)
             +nup.evidence_S_IO/(nup.evidence_S_NO+nup.evidence_S_IO)*nup.fitting_function_S_IO(x), label='Seesaw NO+IO' );
plt.xlabel('$\\Sigma m_\\nu \,\,[\\, {\\rm eV}\\,]$');
plt.ylabel('$\\pi(  \\Sigma m_\\nu) \,\,[\\, {\\rm eV}^{-1}\\,]$');
plt.legend();
plt.savefig(outdir+'1_figure7.pdf')
plt.clf()
