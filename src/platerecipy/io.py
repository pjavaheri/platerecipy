"""
File brief
----------
`io.py`

Module for input/output functions.

This is a part of `platerecipy` package. For license and citation, please
refer to the main repository:
[github.com/pjavaheri/platerecipy](github.com/pjavaheri/platerecipy)

Author(s): 
Pejvak Javaheri; [pejvak.javaheri@mail.utoronto.ca](mailto:pejvak.javaheri@mail.utoronto.ca)
"""


import logging
log = logging.getLogger(__name__)

from . import _INT
from .model import PlateModel
from .grid import SphericalGrid, PartialSphericalGrid, CustomGrid

from scipy.io import savemat
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


import pyvista as pv

def save_as_vtu(
    model               : PlateModel,
    other_fields        = None,
    other_fields_names  = None,
    filename            = 'platerecipy_output.vtu'
):
    """
    Assumes the existance of a mesh object.
    """
    if isinstance(model.grid, SphericalGrid):
        mask = np.linspace(0, model.grid.xs.size-1, model.grid.xs.size).reshape(model.grid.xs.shape)
        mask = (mask >= model.grid.xs.shape[1]-1)&(mask<=(model.grid.xs.shape[0]-1)*model.grid.xs.shape[1])
    else:
        mask = np.ones(model.grid.xs.shape, dtype=np.bool)

    if other_fields is not None:
        for i in range(len(other_fields)):
            model.grid._mesh[
                other_fields_names[i]
            ] = other_fields[i][mask]
    
    model.grid._mesh['plate_IDs']       = model.plate_IDs[mask]
    model.grid._mesh['ID_probs']        = model.ID_probs[mask]
    model.grid._mesh['stacked_field']   = model.stacked_field[mask]
    model.grid._mesh['markers']         = model.markers[mask]

    model.grid._mesh.save(filename)

def save_as_mat(
    model               : PlateModel,
    other_fields        = None,
    other_fields_names  = None,
    filename            = 'platerecipy_output.mat'
) -> None:
    """
    To save the interpolated model as OCTAVE/MATLAB file.

    Parameters
    ----------
    model : PlateModel,
        An input model whose plates are identified.

    other_fields : list, optional,
        A list of additional fields to be included in the output.

    other_fields_names : list, optional,
        A list of field identifiers corresponding to `other_fields`.
    
    filename : str, default='platerecipy_output.mat'
        Output filename.
    """
    
    data = {
        'plate_IDs'     : model.plate_IDs,
        'stacked_field' : model.stacked_field,
        'ID_probs'      : model.ID_probs,
        'markers'       : model.markers,
        'xs'            : model.grid.xs,
        'ys'            : model.grid.ys,
        'zs'            : model.grid.zs
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
    model               : PlateModel,
    other_fields        = None,
    other_fields_names  = None,
    map_to_original     = False,
    filename            = 'platerecipy_output.csv'
):
    """
    To save the interpolated model as a CSV file.

    Parameters
    ----------
    model : PlateModel,
        An input model whose plates are identified.

    other_fields : list, optional,
        A list of additional fields to be included in the output.

    other_fields_names : list, optional,
        A list of field identifiers corresponding to `other_fields`.

    map_to_original : bool, default=False
        Whether to map the fields back to the original input data points.
    
    filename : str, default='platerecipy_output.csv'
        Output filename.
    """
    fields_to_write = [model.plate_IDs, model.stacked_field, model.ID_probs, model.markers]
    fields_names_to_write = ['plate_IDs', 'ID_probs', 'stacked_field', 'markers']

    if other_fields is not None:
        for i in range(len(other_fields)):
            fields_to_write.append(other_fields[i])
            fields_names_to_write.append(other_fields_names[i])

    if map_to_original:
        fields_to_write = [
            model.grid.map_to_original_input(f.ravel(), method='tangent-plane') \
                for f in fields_to_write
        ]
        fields_to_write.append(model.grid.original_xs)
        fields_to_write.append(model.grid.original_ys)
        fields_to_write.append(model.grid.original_zs)
        fields_names_to_write += ["original_xs", "original_ys", "original_zs"]
    else:
        fields_to_write = [f.ravel() for f in fields_to_write]
        fields_to_write.append(model.grid.xs.ravel())
        fields_to_write.append(model.grid.ys.ravel())
        fields_to_write.append(model.grid.zs.ravel())
        fields_names_to_write += ["xs", "ys", "zs"]
    
    fmtstr = ""
    for i in range(len(fields_to_write)):
        if fields_to_write[i].dtype == _INT:
            fmtstr += "{" + str(i) + ":d},"
        else:
            fmtstr += "{" + str(i) + ":.7e},"
    fmtstr = fmtstr[:-1] + '\n'

    with open(filename, 'w') as outfile:
        header = ""
        for var in fields_names_to_write:
            header += f"{var},"
        header = header[:-1] + '\n'
        outfile.write(header)
        for i in range(fields_to_write[0].size):
            outfile.write(
                fmtstr.format(
                    *tuple(arr[i] for arr in fields_to_write)
                )
            )



def save_mollweide_projection(
    model       : PlateModel,
    filename    = 'platerecipy_mollweide_output.png'
) -> None:
    """
    To save the model as a Mollweide projection.

    Parameters
    ----------
    model : PlateModel,
        An input model whose plates are identified.
    
    filename : str, default='platerecipy_mollweide_output.png'
        Output filename.
    """
    
    if not isinstance(model.grid, SphericalGrid):
        raise ValueError('Bad input. A regular spherical grid is required for generating a mollweide projection.')

    with plt.rc_context({
        'mathtext.fontset'  : 'stix', 
        'font.family'       : 'STIXGeneral', 
        'font.size'         : 16
    }):
        fig, axes = plt.subplots(
            3, 1, figsize=(10, 12), subplot_kw={'projection': 'mollweide'}
        )

        # ~~~~~~~~~ for the stacked field ~~~~~~~~~
        ax = axes[0]
        pc = ax.pcolormesh(
            model.grid.phis, np.pi/2 - model.grid.thetas, model.stacked_field, shading='nearest'
        )
        cb = fig.colorbar(pc, ax=ax, orientation="vertical", shrink=0.5)
        cb.set_label("Stacked field")
        ax.set_yticks([-np.pi*5/12, -np.pi/3, -np.pi/4, -np.pi/6, -np.pi/12, 0, np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi*5/12])
        ax.set_xticks([-np.pi*5/6,-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2, np.pi*5/6])


        # ~~~~~~~~ for the plate IDs ~~~~~~~~~
        ax = axes[1]
        bounds = np.linspace(0.5, model.plate_IDs.max()+0.5, model.plate_IDs.max()+1)
        norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pc = ax.pcolormesh(
            model.grid.phis, np.pi/2 - model.grid.thetas, model.plate_IDs, shading='nearest',
            cmap='jet', norm=norm
        )

        if model.plate_IDs.min() == 0:
            # non-conforming regions were identified
            # they will be plotted in black shade
            temp = model.plate_IDs.copy().astype(float)
            temp[model.plate_IDs > 0] = float('NaN')
            ax.pcolormesh(
                model.grid.phis, np.pi/2 - model.grid.thetas, temp, shading='nearest',
                cmap='Grays', vmin=-1, vmax=0
            )
        
        cb = fig.colorbar(pc, ax=ax, orientation="vertical", shrink=0.5)
        cb.set_label("Plate ID")
        cb.ax.set_yticks(
            [1, model.plate_IDs.max()], 
            [1, model.plate_IDs.max()]
        )
        ax.set_yticks([-np.pi*5/12, -np.pi/3, -np.pi/4, -np.pi/6, -np.pi/12, 0, np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi*5/12])
        ax.set_xticks([-np.pi*5/6,-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2, np.pi*5/6])


        # ~~~~~~~~ for ID probs ~~~~~~~~~~~
        ax = axes[2]
        pc = ax.pcolormesh(
            model.grid.phis, np.pi/2 - model.grid.thetas, model.ID_probs, shading='nearest', 
            cmap='coolwarm', vmin=0, vmax=1 
        )
        cb = fig.colorbar(pc, ax=ax, orientation="vertical", shrink=0.5)
        cb.set_label("Segmentation probability")
        ax.set_yticks([-np.pi*5/12, -np.pi/3, -np.pi/4, -np.pi/6, -np.pi/12, 0, np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi*5/12])
        ax.set_xticks([-np.pi*5/6,-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2, np.pi*5/6])
        
        fig.tight_layout()
        fig.savefig(filename, dpi=300)
        plt.close()


def save_six_view_angles(
    model       : PlateModel,
    filename    = 'platerecipy_six_angle_output.png'
):
    """
    To save the model viewed in six angles.

    Parameters
    ----------
    model : PlateModel,
        An input model whose plates are identified.
    
    filename : str, default='platerecipy_six_angle_output.png'
        Output filename.
    """
    if type(model.grid) is PartialSphericalGrid:
        raise ValueError('Bad input. Only a full spherical model can have six angle view output.')
    elif isinstance(model.grid, CustomGrid):
        log.warning("For a six view angle image, the grid must represent the entire spherical surface.")
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    if isinstance(model.grid, SphericalGrid):
        mask = np.linspace(0, model.grid.xs.size-1, model.grid.xs.size).reshape(model.grid.xs.shape)
        mask = (mask >= model.grid.xs.shape[1]-1)&(mask<=(model.grid.xs.shape[0]-1)*model.grid.xs.shape[1])
    else:
        mask = np.ones(model.grid.xs.shape, dtype=np.bool)

    model.grid._mesh['plate_IDs']       = model.plate_IDs[mask]
    model.grid._mesh['ID_probs']        = model.ID_probs[mask]
    model.grid._mesh['stacked_field']   = model.stacked_field[mask]

    mesh = model.grid._mesh

    with plt.rc_context({
        'mathtext.fontset'  : 'stix', 
        'font.family'       : 'STIXGeneral', 
        'font.size'         : 16
    }):
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
                line_width      = 4,
                cone_radius     = 0.5,
                shaft_length    = 1,
                tip_length      = 0.4,
                ambient         = 0.0,
                label_size      = (0.2, 0.2),
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
            pl.window_size = [800, 800]
            pl.add_mesh(mesh, style='surface', line_width=4, cmap='nipy_spectral', scalars='plate_IDs')
            pl.camera_position = camera_positions[i]
            #pl.show_axes()
            pl.add_axes(
                line_width      = 4,
                cone_radius     = 0.5,
                shaft_length    = 1,
                tip_length      = 0.4,
                ambient         = 0.0,
                label_size      = (0.2, 0.2),
            )
            pl.zoom_camera(1.5)
            pl.remove_scalar_bar()
            if i == 5:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes('bottom', size='5%', pad=0.05)
                #fig.colorbar(im, cax=cax, cmap='coolwarm', orientation='horizontal', norm)
                bounds = np.linspace(
                    mesh['plate_IDs'].min(), 
                    mesh['plate_IDs'].max(), 
                    mesh['plate_IDs'].max() - mesh['plate_IDs'].min() + 1
                )
                norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
                c = ax.imshow(pl.screenshot(), cmap='nipy_spectral', norm=norm)
                cbar = fig.colorbar(c, cax=cax, orientation='horizontal')
                cbar.ax.set_xticks([
                    mesh['plate_IDs'].min(), 
                    mesh['plate_IDs'].max()
                ])
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
                line_width      = 4,
                cone_radius     = 0.5,
                shaft_length    = 1,
                tip_length      = 0.4,
                ambient         = 0.0,
                label_size      = (0.2, 0.2),
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
        plt.close()