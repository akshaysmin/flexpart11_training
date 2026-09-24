"""
Module: inverse.functions
Description:
    Utility functions for reading FLEXPART inverse model output (netCDF), creating
    plots of emission sensitivities and time series, manipulating emission
    source definitions, and performing a simple Bayesian inversion to estimate
    emissions from observations.
    This module is intended for diagnostic, visualization, and small-scale
    inversion experiments using sensitivity grids produced by FLEXPART
    (or similar transport models). It assumes specific variable names and
    array shapes in the netCDF files (see read_grid_time_file).
Copyright:
    Copyright (c) 2026, Flexpart Training Contributors.
    Author: Martin Vojta, Michael Blaschek
    All rights reserved.
Functions:
    add_source(emissions: list, *, lat: float, lon: float, val: float, edgecolor="green")
        Add a single emission source to a list of emission descriptors.
        Parameters:
            emissions (list): Mutable list of emission dicts to append to. Each
                emission dict contains keys 'lat', 'lon', 'val', and optional
                'edgecolor'.
            lat (float): Latitude of the emission cell center (degrees).
            lon (float): Longitude of the emission cell center (degrees).
            val (float): Emission value (units used elsewhere in the code;
                stored internally scaled by 1e-12).
            edgecolor (str, optional): Matplotlib edge color for plotting the
                emission rectangle. Defaults to "green".
        Returns:
            list: The same emissions list with the new source appended.
    read_grid_time_file(grid_time_file)
        Read model grid, time and sensitivity data from a netCDF file.
        Parameters:
            grid_time_file (str or path-like): Path to a netCDF file containing
                required variables:
                  - latitude, longitude (1-D arrays)
                  - time (1-D array)
                  - height (1-D array)
                  - spec001_mr (sensitivity array expected with dimensions
                    [nageclass, pointspec, time, height, lat, lon])
                  - RELSTART, RELEND (relative times used to compute release_times)
        Returns:
            tuple: (lat, lon, time, release_times, height, conc, netcdf_dataset)
                - lat (ndarray): latitude coordinates
                - lon (ndarray): longitude coordinates
                - time (ndarray): raw time values from file
                - release_times (list[datetime]): reference times for each release
                - height (ndarray): model height levels
                - conc (ndarray): sensitivity array (as read from 'spec001_mr')
                - f (netCDF4.Dataset): open Dataset object (caller may use/close)
        Notes:
            - The function derives a reference datetime from the filename
              (expects a YYYYMMDDHHMMSS timestamp near the end of the filename).
            - release_times is computed from RELSTART and RELEND variables.
    plot_sensitivity_for_all_releases(grid_time_file, map_coordinates, emissions, colorbar_limits)
        Plot spatial sensitivity maps for all release times in a single figure
        arranged in a grid of subplots.
        Parameters:
            grid_time_file (str): Path to the netCDF file (forwarded to
                read_grid_time_file).
            map_coordinates (tuple): [lon_min, lon_max, lat_min, lat_max] used
                for axis.set_extent.
            emissions (list): List of emission dicts (see add_source). If not
                empty, rectangular patches are added to mark emission grid cells.
            colorbar_limits (tuple): (vmin, vmax) limits used by a logarithmic
                color normalization for the plotted sensitivities.
        Returns:
            matplotlib.figure.Figure: The generated figure.
        Notes:
            - Uses a custom colormap that sets the first color to white to
              emphasize zero/near-zero values.
            - Assumes the sensitivity to plot is conc[0, i, :, 0, :, :] summed
              over the vertical dimension (or as appropriate for the file).
    create_rectangle(lat, lon, emissions, edgecolor="green")
        Create matplotlib Rectangle patches centered on grid cells corresponding
        to emission locations.
        Parameters:
            lat (ndarray): 1-D latitude coordinates.
            lon (ndarray): 1-D longitude coordinates.
            emissions (list): List of emission dicts (each with 'lat' and 'lon').
            edgecolor (str, optional): Default edge color for rectangles.
        Returns:
            list[matplotlib.patches.Rectangle]: List of rectangle patch objects
            using PlateCarree transform ready to be added to GeoAxes.
    calculate_timeseries(lat, lon, height, conc, emissions)
        Compute a synthetic concentration time series from sensitivity grids and
        a set of point emissions.
        Parameters:
            lat, lon (ndarray): Grid coordinate arrays.
            height (ndarray): Model heights (expects height[0] to be scalar or
                contain a normalization factor).
            conc (ndarray): Sensitivity array with shape consistent with
                read_grid_time_file output.
            emissions (list): List of emission dicts with keys 'lat','lon','val'.
        Returns:
            ndarray: 1-D timeseries (same length as number of releases / times)
                expressed in ppt (part-per-trillion) given the internal unit
                conventions in this code:
                  - sensitivity (sm^3/kg) * flux (ng/m^2/s) / height => ng/m^3/s
                  - multiplied by 1e12 to yield ppt in the code's convention.
        Notes:
            - The function places the emission value at the nearest grid cell
              index for each source and multiplies element-wise with the summed
              sensitivities to produce contributions to the timeseries.
    plot_timeseries(release_times, timeseries, label=None, ax=None)
        Simple helper to plot a timeseries of concentrations vs release_times.
        Parameters:
            release_times (list[datetime] or ndarray): Time points for x-axis.
            timeseries (ndarray): Concentration values to plot.
            label (str, optional): Line label for legend.
            ax (matplotlib.axes.Axes, optional): Axes to draw into; if None a
                new figure/axes pair is created.
        Returns:
            matplotlib.axes.Axes: Axes containing the plotted series.
    perturb_emissions(pertubations, emissions)
        Build a perturbed (a priori) emissions vector from fractional perturbations.
        Parameters:
            pertubations (iterable): Relative perturbation factors (e.g. [+0.1, -0.2]).
            emissions (list): List of emission dicts (with 'val' keys).
        Returns:
            ndarray: a_priori_emissions (1-D array) computed as emissions[i]['val'] *
            (1 + pertubations[i]).
        Notes:
            - The function name and parameter name 'pertubations' follow the
              calling code; ensure correct spelling when invoking.
    concentration_and_emissions_before_inversion(lat, lon, height, release_times, conc, e, xp)
        Visual comparison of true vs prior emissions and the corresponding
        concentration time series before any inversion is applied.
        Parameters:
            lat, lon, height, release_times, conc: Geometry and sensitivity arrays
                as used elsewhere.
            e (list): True emission dicts.
            xp (iterable): Prior emission values (same length as e).
        Returns:
            matplotlib.figure.Figure: Figure containing two subplots:
                - concentration time series for prior and true emissions
                - bar plot comparing prior and true emission values (ng/m^2/s)
    stats(y_true, y_pred)
        Print simple error statistics comparing two 1-D arrays:
            - MSE, RMSE, MAE, MAPE (percentage)
        Parameters:
            y_true (ndarray): Reference/true values.
            y_pred (ndarray): Predicted/estimated values.
        Returns:
            None (prints results).
        Notes:
            - MAPE divides by y_true; ensure no zeros in y_true when calling.
    concentration_and_emissions_after_inversion(lat, lon, height, release_times, conc, e, xp_prior, xp_post)
        Visualize and compare true, prior and posterior emissions and their
        concentrations after an inversion.
        Parameters:
            lat, lon, height, release_times, conc: Geometry and sensitivity arrays.
            e (list): True emission dicts.
            xp_prior (iterable): Prior emission values.
            xp_post (iterable): Posterior emission values (inferred).
        Returns:
            matplotlib.figure.Figure: Figure showing:
                - concentration time series for prior, posterior and true emissions
                - bar plot comparing prior, posterior and true emission values
        Side-effects:
            - Prints statistics (MSE, RMSE, MAE, MAPE) comparing prior/posterior
              time series to true time series.
    calculate_Transport_matrix_H(lat, lon, conc, height, e)
        Construct the transport/sensitivity matrix H (and its transpose Ht)
        mapping emissions (grid cells/sources) to observations (release times).
        Parameters:
            lat, lon (ndarray): Coordinate arrays.
            conc (ndarray): Sensitivity array as returned by read_grid_time_file.
            height (ndarray): Height array (used for normalization).
            e (list): List of emission dicts (with lat/lon locations identifying
                grid cells that correspond to columns of H).
        Returns:
            tuple: (H, Ht)
                - H (ndarray): Observation-by-source matrix (time x n_sources)
                - Ht (ndarray): Transpose of H (n_sources x time)
        Notes:
            - The function averages/sums over the model vertical levels using
              height[0] and selects grid cells corresponding to each source.
            - The implementation constructs H as time-major and then returns
              both H and its transpose.
    inversion(xp, emission_error, observation_error, timeseries, H, Ht)
        Perform a simple linear Bayesian (Gauss-Markov) inversion to estimate
        posterior emissions from observations.
        Parameters:
            xp (ndarray): Prior emission values (state vector).
            emission_error (ndarray or iterable): Diagonal elements (relative
                variances or factors) used to form the background covariance B;
                B = diag(emission_error * xp). Units must be consistent with xp.
            observation_error (float or scalar-like): Measurement error value
                used to build R as observation_error * I (applied to each time).
            timeseries (ndarray): Observations (e.g. concentration timeseries)
                in the same unit convention expected by H and xp.
            H (ndarray): Observation-by-source transport matrix.
            Ht (ndarray): Transpose of H (source-by-observation).
        Returns:
            ndarray: Posterior emission estimate (same shape as xp).
        Implementation notes:
            - The function applies small-unit scaling: it converts observation
              and error quantities by 1e-12 before inversion to match the
              internal unit conventions in this module.
            - The posterior update uses the standard formula:
                x_post = xp + B H^T (H B H^T + R)^{-1} (y - H xp)
            - R is constructed as diag(observation_error repeated) * 1e-12.
            - B is constructed as diag(emission_error * xp).
            - Users should ensure H, xp, and error specifications are consistent
              and well-conditioned; no explicit regularization other than B
              and R is applied.
General Notes and Assumptions:
    - Many functions assume 0-based indexing and that lat/lon arrays are
      regularly spaced so that lat[1] - lat[0] and lon[1] - lon[0] yield grid
      cell extents.
    - Units and scaling: the module contains scaling factors (notably 1e-12
      and 1e12) that are chosen to yield convenient plotting units (ppt,
      ng/m^2/s, etc.). Verify consistency for your datasets.
    - The netCDF reading logic expects particular variable names; adapt
      read_grid_time_file if using different file conventions.
    - The plotting helpers rely on cartopy and matplotlib; make sure a proper
      PROJ/data environment is available for cartopy features to render.
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors

#######################################################################################################################################################

def add_source(emissions: list, *, lat: float, lon: float, val: float, edgecolor="green"):
    emissions.append({"lat": lat, "lon": lon, "val": val * 1e-12, "edgecolor": edgecolor})
    return emissions


########################################################################################################################################################

def read_grid_time_file(grid_time_file):
    from datetime import datetime, timedelta

    f = Dataset(grid_time_file)
    lat = f.variables["latitude"][:]
    lon = f.variables["longitude"][:]
    time = f.variables["time"][:]
    height = f.variables["height"][:]
    conc = f.variables["spec001_mr"][
        :, :, :, :, :, :
    ]  # nageclass, pointspec, time, height, latitude, longitude

    t = grid_time_file[-17:-3]
    ref_time = datetime.strptime(t, "%Y%m%d%H%M%S")
    relative_seconds = (f.variables["RELSTART"][:] + f.variables["RELEND"][:]) / 2
    release_times = [ref_time + timedelta(seconds=int(sec)) for sec in relative_seconds]

    return lat, lon, time, release_times, height, conc, f


###############################################################################################################################################################


def plot_sensitivity_for_all_releases(
    grid_time_file, map_coordinates, emissions, colorbar_limits
):

    lat, lon, time, release_times, height, conc, f = read_grid_time_file(grid_time_file)

    colors = plt.cm.RdPu(np.linspace(0, 1, 256))
    colors[0] = [1, 1, 1, 1]
    cmap = ListedColormap(colors)

    releases = len(conc[0])
    ncols = int(np.ceil(np.sqrt(releases)))
    nrows = int(np.ceil(releases / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(4 * ncols, 4 * nrows),
    )

    axes = np.array(axes).flatten()

    for i in range(releases):
        ax = axes[i]
        # sum over all time steps and/or heights too? (could be added here)
        # ageclass, pointspec, time, height, lat, lon
        data = np.sum(conc[0, i, :, 0, :, :], axis=0)
        im = ax.pcolormesh(
            lon,
            lat,
            data,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            norm=mcolors.LogNorm(vmin=colorbar_limits[0], vmax=colorbar_limits[1]),
        )
    
        ax.set_title(f"t = {release_times[i]}")
        ax.set_extent(map_coordinates)
        ax.add_feature(cartopy.feature.COASTLINE.with_scale("110m"), alpha=0.5)
        ax.add_feature(cartopy.feature.BORDERS.with_scale("110m"), alpha=0.5)
        ax.gridlines(draw_labels=False)
        if len(emissions) > 0:
            rects = create_rectangle(lat, lon, emissions)
            for rect in rects:
                ax.add_patch(rect)

        gl = ax.gridlines(draw_labels=True, linestyle=":", alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False

        # Only label left column and bottom row
        row = i // ncols
        col = i % ncols
        gl.left_labels = col == 0
        gl.bottom_labels = row == nrows - 1

    for j in range(releases, len(axes)):
        fig.delaxes(axes[j])

    cbar_ax = fig.add_axes([1.01, 0.1, 0.03, 0.80])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar.set_label("Emission sensitivity [sm$^3$/kg]", fontsize=20)
    cbar.ax.tick_params(labelsize=16)

    plt.show()
    return fig


##############################################################################################################


def create_rectangle(lat, lon, emissions, edgecolor="green"):
    from matplotlib.patches import Rectangle

    rects = []
    for i in range(len(emissions)):
        target_lat = emissions[i]["lat"]
        target_lon = emissions[i]["lon"]
        ind_lat = np.argmin(np.abs(lat - target_lat))
        ind_lon = np.argmin(np.abs(lon - target_lon))
        latbox = lat[ind_lat]
        lonbox = lon[ind_lon]
        dlat = lat[1] - lat[0]
        dlon = lon[1] - lon[0]
        ll_lat = latbox - dlat / 2
        ll_lon = lonbox - dlon / 2
        rect = Rectangle(
            (ll_lon, ll_lat),
            dlon,
            dlat,
            linewidth=3,
            edgecolor=emissions[i].get("edgecolor", edgecolor),
            facecolor="none",
            transform=ccrs.PlateCarree(),
        )
        rects.append(rect)
    return rects


#######################################################################################


def calculate_timeseries(lat, lon, height, conc, emissions):

    em = np.zeros((len(lat), len(lon)))
    for i in range(len(emissions)):
        target_lat = emissions[i]["lat"]
        target_lon = emissions[i]["lon"]
        em_val = emissions[i]["val"]
        ind_lat = np.argmin(np.abs(lat - target_lat))
        ind_lon = np.argmin(np.abs(lon - target_lon))
        em[ind_lat, ind_lon] = em_val

    timeseries = np.zeros(len(conc[0]))
    for i in range(len(timeseries)):
        # Compute concentration weighted by emissions
        # This requires ind_receptor to be 2 to get sm3/kg as units
        # s m3 / kg * ng / m2 / s * 1e12 = ppt (part per trillion)
        c = np.sum(conc[0, i, :, 0, :, :] / height[0], axis=0) * em
        timeseries[i] = np.sum(c) * 10**12
    return timeseries


#############################################################################################


def plot_timeseries(release_times, timeseries, label=None, ax=None):

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))
    ax.plot(release_times, timeseries, label=label)
    ax.set_title("concentration time series")
    ax.set_xlabel("Release Time")
    ax.set_ylabel("Concentration [ppt]")
    ax.tick_params(axis="x", labelrotation=45)
    return ax


#############################################################################################


def perturb_emissions(pertubations, emissions):
    a_priori_emissions = np.array(
        [
            emissions[i]["val"] + pertubations[i] * emissions[i]["val"]
            for i in range(len(emissions))
        ]
    )
    return a_priori_emissions


###############################################################################################


def concentration_and_emissions_before_inversion(
    lat, lon, height, release_times, conc, e, xp
):

    import numpy as np
    import matplotlib.pyplot as plt

    # Compute time series for true and perturbed emissions
    timeseries = calculate_timeseries(lat, lon, height, conc, e)

    ep = [src.copy() for src in e]
    for i in range(len(ep)):
        ep[i]["val"] = xp[i]

    timeseries_perturbed = calculate_timeseries(lat, lon, height, conc, ep)

    true_emissions = np.array([e[i]["val"] for i in range(len(e))])
    perturbed_emissions = xp

    a = perturbed_emissions * 1e12
    b = true_emissions * 1e12

    indices = np.arange(len(a))
    width = 0.35

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(release_times, timeseries_perturbed, label="prior emissions", color="blue")
    ax1.plot(release_times, timeseries, label="true emissions", color="orange")
    ax1.set_title("Concentration Time Series")
    ax1.set_xlabel("Release Time")
    ax1.set_ylabel("Concentration [ppt]")
    ax1.tick_params(axis="x", labelrotation=45)
    ax1.legend()
    ax1.grid()

    ax2.bar(indices - width / 2, a, width=width, label="prior emissions", color="blue")
    ax2.bar(indices + width / 2, b, width=width, label="true emissions", color="orange")
    ax2.set_xlabel("emission grid cell")
    ax2.set_ylabel("emissions [ng/m²/s]")
    ax2.set_title("Emission Values by grid cell")
    ax2.set_xticks(indices)
    ax2.legend()

    plt.tight_layout()
    plt.show()

    return fig



def stats(y_true, y_pred):
    import numpy as np
    # Mean Squared Error (MSE)
    mse = np.mean((y_true - y_pred) ** 2)

    # Root Mean Squared Error (RMSE)
    rmse = np.sqrt(mse)

    # Mean Absolute Error (MAE)
    mae = np.mean(np.abs(y_true - y_pred))

    # Mean Absolute Percentage Error (MAPE)
    mape = np.mean(np.abs((y_true - y_pred) / y_true)) * 100

    print(f"MSE:  {mse:.2f}")
    print(f"RMSE: {rmse:.2f}")
    print(f"MAE:  {mae:.2f}")
    print(f"MAPE: {mape:.2f}%")

#####################################################################################################


def concentration_and_emissions_after_inversion(
    lat, lon, height, release_times, conc, e, xp_prior, xp_post
):
    import numpy as np
    import matplotlib.pyplot as plt

    # Compute time series for true emissions
    timeseries_true = calculate_timeseries(lat, lon, height, conc, e)

    e_prior = [src.copy() for src in e]
    e_post = [src.copy() for src in e]

    for i in range(len(e)):
        e_prior[i]["val"] = xp_prior[i]
        e_post[i]["val"] = xp_post[i]

    timeseries_prior = calculate_timeseries(lat, lon, height, conc, e_prior)
    timeseries_post = calculate_timeseries(lat, lon, height, conc, e_post)

    print("Statistics for Prior vs True Emissions:")
    stats(timeseries_true, timeseries_prior)
    print("Statistics for Posterior vs True Emissions:")
    stats(timeseries_true, timeseries_post)

    true_emissions = np.array([e[i]["val"] for i in range(len(e))])
    prior_emissions = xp_prior
    posterior_emissions = xp_post

    scale = 1e12
    a = prior_emissions * scale
    b = posterior_emissions * scale
    c = true_emissions * scale

    indices = np.arange(len(a))
    width = 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

    ax1.plot(release_times, timeseries_prior, label="prior emissions", color="blue")
    ax1.plot(release_times, timeseries_post, label="posterior emissions", color="green")
    ax1.plot(release_times, timeseries_true, label="true emissions", color="orange")

    ax1.set_title("Concentration Time Series")
    ax1.set_xlabel("Release Time")
    ax1.set_ylabel("Concentration [ppt]")
    ax1.tick_params(axis="x", labelrotation=45)
    ax1.legend()
    ax1.grid()

    ax2.bar(indices - width, a, width=width, label="prior emissions", color="blue")
    ax2.bar(indices, b, width=width, label="posterior emissions", color="green")
    ax2.bar(indices + width, c, width=width, label="true emissions", color="orange")

    ax2.set_xlabel("Emission Grid Cell")
    ax2.set_ylabel("Emissions [ng/m²/s]")
    ax2.set_title("Emission Values by Grid Cell")
    ax2.set_xticks(indices)
    ax2.legend()

    plt.tight_layout()
    plt.show()

    return fig


#################################################################################################################


def calculate_Transport_matrix_H(lat, lon, conc, height, e):
    H = []
    for i in range(len(e)):
        target_lat = e[i]["lat"]
        target_lon = e[i]["lon"]
        ind_lat = np.argmin(np.abs(lat - target_lat))
        ind_lon = np.argmin(np.abs(lon - target_lon))
        H.append(
            np.sum(
                np.array(
                    [conc[0, i, :, 0, ind_lat, ind_lon] for i in range(len(conc[0]))]
                ),
                axis=1,
            )
        )
    H = H / height[0]
    H = H.T
    Ht = np.array(H).T
    return H, Ht


##########################################################################################################


def inversion(xp, emission_error, observation_error, timeseries, H, Ht):
    R = np.diag(np.repeat(observation_error, len(timeseries))) * 1e-12
    # Ri = np.linalg.inv(R)
    B = np.diag(emission_error * xp)
    # Bi = np.linalg.inv(B)
    y = np.array(timeseries) * 1e-12
    # Compute the posterior emissions using the inversion formula
    x = xp + B @ Ht @ (np.linalg.inv(H @ B @ Ht + R)) @ (y - H @ xp)
    return x
