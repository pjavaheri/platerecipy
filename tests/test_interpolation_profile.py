"""
@file test_interpolation_profile.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Interpolation method profiles.
"""


import numpy as np
import matplotlib.pyplot as plt
from platerecipy.grid import SphericalGrid, PartialSphericalGrid
from scipy.special import sph_harm_y

def get_interp_err(func, method, N, res, R=6400e3, ret=False):
    rng         = np.random.default_rng(103)
    data_thetas = np.pi*rng.random(N)
    data_phis   = 2*np.pi*rng.random(N) - np.pi

    data_xs = R*np.sin(data_thetas)*np.cos(data_phis)
    data_ys = R*np.sin(data_thetas)*np.sin(data_phis)
    data_zs = R*np.cos(data_thetas)
    grid = SphericalGrid(
        data_xs, data_ys, data_zs, 
        theta_res=res, phi_res=2*res,
    )
    data_fld   = func(data_thetas, data_phis)*1e-10
    expect_fld = func(grid.thetas, grid.phis)*1e-10
    interp_fld = grid.interpolate_field(data_fld, method=method, near_thresh=0)
    if ret:
        return expect_fld
    dA = np.sin(grid.thetas)
    expect_fld_avg = np.sqrt(
        np.sum(dA*(expect_fld)**2) / np.sum(dA)
    )
    return np.sqrt(
        np.sum(dA*(interp_fld-expect_fld)**2) / np.sum(dA)
    ) / expect_fld_avg
    
fig, axes = plt.subplots(3, 3, figsize=(12, 10))

func = lambda thetas, phis: sph_harm_y(0, 0, thetas, phis).real
ax = axes[0][0]
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=0}^{m=0}; n_{\theta}=2^8$')

Ns = np.array([100, 1000, 10_000, 100_000])
res = 256
method = 'lat-lon'
[
    ax.plot(
        Ns,
        [
            get_interp_err(func, method, N, res) \
                for N in Ns
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel('$N$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[0][1]
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=0}^{m=0}; N=10^4$')
N = 10_000
ress = [16, 32, 64, 128, 256, 512, 1024]
method = 'lat-lon'
[
    ax.plot(
        ress,
        [
            get_interp_err(func, method, N, res) \
                for res in ress
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel(r'$n_{\theta}=n_{\phi}/2$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[0][2]
ax.set_title("Expected")
fig.colorbar(ax.imshow(get_interp_err(func, method, 10_000, 256, ret=True)), orientation='horizontal')
ax.set_axis_off()

func = lambda thetas, phis: sph_harm_y(2, 1, thetas, phis).real
ax = axes[1][0]
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=2}^{m=1}; n_{\theta}=2^8$')

Ns = np.array([100, 1000, 10_000, 100_000])
res = 256
method = 'lat-lon'
[
    ax.plot(
        Ns,
        [
            get_interp_err(func, method, N, res) \
                for N in Ns
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel(r'$N$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')


ax = axes[1][1]
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=0}^{m=0}; N=10^4$')
N = 10_000
ress = [16, 32, 64, 128, 256, 512, 1024]
method = 'lat-lon'
[
    ax.plot(
        ress,
        [
            get_interp_err(func, method, N, res) \
                for res in ress
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel(r'$n_{\theta}=n_{\phi}/2$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[1][2]
ax.set_title("Expected")
fig.colorbar(ax.pcolormesh(get_interp_err(func, method, 10_000, 256, ret=True), shading='neares'), orientation='horizontal')
ax.set_axis_off()

func = lambda thetas, phis: sph_harm_y(12, 5, thetas, phis).real
ax = axes[2][0]
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=12}^{m=5}; n_{\theta}=2^8$')
Ns = np.array([100, 1000, 10_000, 100_000])
res = 256
method = 'lat-lon'
[
    ax.plot(
        Ns,
        [
            get_interp_err(func, method, N, res) \
                for N in Ns
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel('$N$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[2][1]
N = 10_000
ax.set_title(r'$f = 10^{-10} \mathcal{Y}_{\ell=12}^{m=5}; N=10^4$')
ress = [16, 32, 64, 128, 256, 512]
method = 'lat-lon'
[
    ax.plot(
        ress,
        [
            get_interp_err(func, method, N, res) \
                for res in ress
        ],
        '-o',
        label="lat-lon (scipy's griddata)" if method=='lat-lon' else method
    ) for method in ['lat-lon', 'triangle', 'tangent-plane']
]
ax.set_xlabel(r'$n_{\theta}=n_{\phi}/2$')
ax.set_ylabel("RMSE/RMS")
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[2][2]
ax.set_title("Expected")
fig.colorbar(ax.imshow(get_interp_err(func, method, 10_000, 256, ret=True)), orientation='horizontal')
ax.set_axis_off()



fig.tight_layout()
fig.savefig("test_interpolation_profile.png", dpi=300)