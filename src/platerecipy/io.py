"""
@file io.py
@author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
@brief Module for input/output functions.
"""

from platerecipy.model import PlateModel
from platerecipy.grid import convert_grid_to_mesh
from scipy.io import savemat
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt

def save_as_mat(
    model               : PlateModel,
    other_fields        = None,
    other_fields_names  = None,
    filename            = 'platerecipy_output.mat'
):
    data = {
        'plate_IDs'     : model.plate_IDs,
        'stacked_field' : model.stacked_field
    }

    if other_fields is not None:
        for i in range(len(other_fields)):
            data[other_fields_names[i]] = other_fields[i]
    
    savemat(
        file_name           = filename,
        mdict               = data,
        appendmat           = True,
        format              = '5',
        long_field_names    = False,
        do_compression      = False,
        oned_as             = 'row',
    )

def save_as_csv(
    model       : PlateModel,
    other_fields = None,
    other_fields_names = None,
    filename    = 'platerecipy_output.csv'
):
    data = {
        'plate_IDs'     : model.plate_IDs.ravel(),
        'stacked_field' : model.stacked_field.ravel()
    }
    if other_fields is not None:
        for i in range(len(other_fields)):
            data[other_fields_names[i]] = other_fields[i].ravel()
    pd.DataFrame(data).to_csv(filename)

def save_as_vtp(
    model       : PlateModel,
    other_fields        = None,
    other_fields_names  = None,
    filename            = 'platerecipy_output.vtp',
    return_mesh         = False
):
    if other_fields is not None:
        other_fields.append(model.plate_IDs)
        other_fields.append(model.ID_probs)
        other_fields.append(model.stacked_field)
        other_fields_names.append('plate_IDs')
        other_fields_names.append('ID_probs')
        other_fields_names.append('stacked_field')
        mesh = convert_grid_to_mesh(
            gridded_fields  = other_fields,
            field_names     = other_fields_names,
            radius          = model.R
        )
        mesh.save(filename)
    else:
        mesh = convert_grid_to_mesh(
            gridded_fields  = [model.plate_IDs, model.ID_probs, model.stacked_field],
            field_names     = ['plate_IDs', 'ID_probs', 'stacked_field'],
            radius          = model.R
        )
        mesh.save(filename)
    if return_mesh:
        return mesh
    
def save_six_view_angles(
    model               : PlateModel,
    filename            = 'platerecipy_six_angle_output.png'
):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib as mpl
    import numpy as np
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    mpl.rcParams['font.size'] = 16

    mesh = convert_grid_to_mesh(
        gridded_fields  = [model.plate_IDs, model.ID_probs, model.stacked_field],
        field_names     = ['plate_IDs', 'ID_probs', 'stacked_field'],
        radius          = 1.
    )

    fig, axes = plt.subplots(6, 3, figsize=(9, 18))

    camera_positions = [
        "xy", "xz", "yz", "yx", "zx", "zy", "iso"
    ]

    axes[0][0].set_title("Stacked field")
    axes[0][1].set_title("Plate ID")
    axes[0][2].set_title("Segmentation probability")

    for i in range(6):
        ax = axes[i][0]
        pl = pv.Plotter(off_screen=True)
        pl.window_size = [800,800]
        pl.add_mesh(mesh, style='surface', line_width=4, cmap='viridis', scalars='stacked_field')
        pl.camera_position = camera_positions[i]
        #pl.show_axes()
        pl.add_axes(
            line_width=4,
            cone_radius=0.5,
            shaft_length=1,
            tip_length=0.4,
            ambient=0.0,
            label_size=(0.2, 0.2),
        )
        pl.zoom_camera(1.5)
        pl.remove_scalar_bar()
        if i == 5:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            #fig.colorbar(im, cax=cax, cmap='coolwarm', orientation='horizontal', norm)
            bounds = np.linspace(0, 1, 10)
            norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            c = ax.imshow(pl.screenshot(), cmap='viridis', vmin=0, vmax=1)
            cbar = fig.colorbar(c, cax=cax, orientation='horizontal')
            cbar.ax.set_xticks([0, 1])
        else:
            im = ax.imshow(pl.screenshot())
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        pl.close()
        ax.set_ylabel(f"${camera_positions[i]}$-plane view")
        del pl


    for i in range(6):
        ax = axes[i][1]
        pl = pv.Plotter(off_screen=True)
        pl.window_size = [800,800]
        pl.add_mesh(mesh, style='surface', line_width=4, cmap='nipy_spectral', scalars='plate_IDs')
        pl.camera_position = camera_positions[i]
        #pl.show_axes()
        pl.add_axes(
            line_width=4,
            cone_radius=0.5,
            shaft_length=1,
            tip_length=0.4,
            ambient=0.0,
            label_size=(0.2, 0.2),
        )
        pl.zoom_camera(1.5)
        pl.remove_scalar_bar()
        if i == 5:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            #fig.colorbar(im, cax=cax, cmap='coolwarm', orientation='horizontal', norm)
            bounds = np.linspace(1, mesh['plate_IDs'].max(), mesh['plate_IDs'].max())
            norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            c = ax.imshow(pl.screenshot(), cmap='nipy_spectral', norm=norm)
            cbar = fig.colorbar(c, cax=cax, orientation='horizontal')
            cbar.ax.set_xticks([1, mesh['plate_IDs'].max()])
        else:
            im = ax.imshow(pl.screenshot())
        ax.set_axis_off()
        pl.close()
        del pl


    for i in range(6):
        ax = axes[i][2]
        pl = pv.Plotter(off_screen=True)
        pl.window_size = [800,800]
        pl.add_mesh(mesh, style='surface', line_width=4, cmap='coolwarm', scalars='ID_probs')
        pl.camera_position = camera_positions[i]
        #pl.show_axes()
        pl.add_axes(
            line_width=4,
            cone_radius=0.5,
            shaft_length=1,
            tip_length=0.4,
            ambient=0.0,
            label_size=(0.2, 0.2),
        )
        pl.zoom_camera(1.5)
        pl.remove_scalar_bar()
        if i == 5:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            #fig.colorbar(im, cax=cax, cmap='coolwarm', orientation='horizontal', norm)
            bounds = np.linspace(0, 1, 10)
            norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            c = ax.imshow(pl.screenshot(), cmap='coolwarm', vmin=0, vmax=1)
            cbar = fig.colorbar(c, cax=cax, orientation='horizontal')
            cbar.ax.set_xticks([0, 1])
        else:
            im = ax.imshow(pl.screenshot())
            
        ax.set_axis_off()
        pl.close()
        del pl

    fig.tight_layout()
    fig.savefig(filename, dpi=500)