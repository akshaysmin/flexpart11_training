from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

# import cartopy.feature as cfeature
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
        data = c = np.sum(conc[0, i, :, 0, :, :], axis=0)
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
