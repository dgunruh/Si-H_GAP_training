# general imports
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as cp


from quippy.potential import Potential
from quippy.descriptors import Descriptor

# ase imports
import ase.io
from ase import Atoms, Atom
from ase import units
from ase.build import molecule

import sys


def rms_dict(x_ref, x_pred):
    """ Takes two datasets of the same shape and returns a dictionary containing RMS error data"""

    x_ref = np.array(x_ref)
    x_pred = np.array(x_pred)

    if np.shape(x_pred) != np.shape(x_ref):
        raise ValueError("WARNING: not matching shapes in rms")

    error_2 = (x_ref - x_pred) ** 2

    average = np.sqrt(np.average(error_2))
    std_ = np.sqrt(np.sqrt(np.var(error_2)))

    return {"rmse": average, "std": std_}


def energy_plot(in_file, out_file, title="Energies", plotFolder=""):
    """ Plots the distribution of energy per atom on the output vs the input"""
    # read files
    in_atoms = ase.io.read(in_file, ":")
    ener_in = [at.info["dft_energy"] / len(at) for at in in_atoms]

    ener_in_array = np.array(ener_in)
    ener_in_array = np.delete(ener_in_array, ener_in_array.argmax())

    out_atoms = ase.io.read(out_file, ":")
    # in_atoms = ase.io.extxyz.read_xyz(in_file)
    # out_atoms = ase.io.extxyz.read_xyz(out_file)

    # list energies
    # DFT energy
    # ener_in = [at.get_potential_energy() / len(at.get_chemical_symbols()) for at in in_atoms]
    ener_out = [at.info["energy"] / len(at) for at in out_atoms]
    ener_out_array = np.array(ener_out)
    ener_out_array = np.delete(ener_out_array, ener_out_array.argmax())

    # print(ener_in)
    # print(ener_out)

    # GAP energy
    # gap2b = Potential(param_filename='GAP.xml')
    # ener_out = []
    # for at in out_atoms:
    #    at.set_calculator(gap2b)
    #    ener_out.append(at.get_potential_energy())

    # out_atoms.set_calculator(gap2b)
    # ener_out = [at.get_potential_energy() / len(at.get_chemical_symbols()) for at in out_atoms]

    # scatter plot of the data
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(7, 5)
    ax.scatter(ener_in_array, ener_out_array)
    # get the appropriate limits for the plot
    for_limits = np.concatenate((ener_in_array, ener_out_array))
    elim = (np.amin(for_limits) - 0.5, np.amax(for_limits) + 0.5)
    ax.set_xlim(elim)
    ax.set_ylim(elim)
    # add line of slope 1 for refrence
    ax.plot(elim, elim, c="k")
    # set labels
    ax.set_ylabel("GAP energy [eV]")
    ax.set_xlabel("DFT energy [eV]")
    # set title
    ax.set_title(title)
    # add text about RMSE
    _rms = rms_dict(ener_in_array, ener_out_array)
    rmse_text = (
        "RMSE:\n" + str(np.round(_rms["rmse"], 3)) + " eV/atom"
    )  # + ' +- ' + str(np.round(_rms['std'], 3)) + ' eV/atom'
    ax.text(
        0.9,
        0.1,
        rmse_text,
        transform=ax.transAxes,
        fontsize="large",
        horizontalalignment="right",
        verticalalignment="bottom",
    )
    plt.savefig(plotFolder + "Energies_SOAP.png", format="png", dpi=300)
    plt.show()


def stress_plot(in_file, out_file, title="Stresses:", plotFolder=""):
    """ Plots the distribution of stress on the output vs the input"""
    # read files
    in_atoms = ase.io.read(in_file, ":")
    stress_in = [at.info["dft_virial"] for at in in_atoms]
    stress_in_array = np.asarray(stress_in)
    # print(stress_in_array)
    # mag_stress_in = [np.linalg.norm(i) for i in stress_in]

    out_atoms = ase.io.read(out_file, ":")
    stress_out = [np.ravel(at.info["virial"]) for at in out_atoms]
    stress_out_array = np.asarray(stress_out)
    # print(stress_out_array)
    # mag_stress_out = [np.linalg.norm(i) for i in stress_out]

    # print(ener_in)
    # print(ener_out)

    # GAP energy
    # gap2b = Potential(param_filename='GAP.xml')
    # ener_out = []
    # for at in out_atoms:
    #    at.set_calculator(gap2b)
    #    ener_out.append(at.get_potential_energy())

    # out_atoms.set_calculator(gap2b)
    # ener_out = [at.get_potential_energy() / len(at.get_chemical_symbols()) for at in out_atoms]

    # scatter plot of the data
    fig, ax_list = plt.subplots(3, 3, gridspec_kw={"hspace": 0.5, "wspace": 0.5})
    fig.set_size_inches(21, 15)
    ax_list = ax_list.flat[:]

    components = [
        r"$\sigma_{xx}$",
        r"$\sigma_{xy}$",
        r"$\sigma_{xz}$",
        r"$\sigma_{yx}$",
        r"$\sigma_{yy}$",
        r"$\sigma_{yz}$",
        r"$\sigma_{zx}$",
        r"$\sigma_{zy}$",
        r"$\sigma_{zz}$",
    ]

    for n, ax in enumerate(ax_list):
        mag_stress_in = stress_in_array[:, n]
        mag_stress_out = stress_out_array[:, n]

        # print(mag_stress_in)
        # print(mag_stress_out)

        ax.scatter(mag_stress_in, mag_stress_out)
        # get the appropriate limits for the plot
        for_limits = np.concatenate((mag_stress_in, mag_stress_out))
        stresslim = (for_limits.min() - 0.05, for_limits.max() + 0.05)
        ax.set_xlim(stresslim)
        ax.set_ylim(stresslim)
        # add line of slope 1 for refrence
        ax.plot(stresslim, stresslim, c="k")
        # set labels
        ax.set_ylabel("GAP virial stress [eV]")
        ax.set_xlabel("DFT virial stress [eV]")
        # set title
        ax.set_title("Stress tensor component: " + components[n])
        # add text about RMSE
        _rms = rms_dict(mag_stress_in, mag_stress_out)
        rmse_text = (
            "RMSE:\n" + str(np.round(_rms["rmse"], 3)) + " eV"
        )  # + ' +- ' + str(np.round(_rms['std'], 3)) + ' eV'
        ax.text(
            0.9,
            0.1,
            rmse_text,
            transform=ax.transAxes,
            fontsize="large",
            horizontalalignment="right",
            verticalalignment="bottom",
        )

    plt.savefig(plotFolder + "Stresses_SOAP.png", format="png", dpi=300)
    plt.show()


def force_plot(in_file, out_file, symbol="HSi", title="Forces", plotFolder=""):
    """Plots the distribution of force components per atom on the output vs the input
    only plots for the given atom type(s)"""

    in_atoms = ase.io.read(in_file, ":")
    out_atoms = ase.io.read(out_file, ":")

    # extract data for only one species
    in_force, out_force = [], []
    for at_in, at_out in zip(in_atoms, out_atoms):
        # get the symbols
        sym_all = at_in.get_chemical_symbols()
        # add force for each atom
        for j, sym in enumerate(sym_all):
            if sym in symbol:
                # in_force.append(at_in.get_forces()[j])
                in_force.append(at_in.arrays["dft_force"][j])
                # out_force.append(at_out.get_forces()[j]) \
                out_force.append(
                    at_out.arrays["force"][j]
                )  # because QUIP and ASE use different names
    # convert to np arrays, much easier to work with
    # in_force = np.array(in_force)
    # out_force = np.array(out_force)
    # scatter plot of the data

    in_force_array = np.asarray(in_force)
    out_force_array = np.asarray(out_force)

    components = [r"$\hat{x}$", r"$\hat{y}$", r"$\hat{z}$"]
    fig, ax_list = plt.subplots(3, 1, gridspec_kw={"hspace": 0.5})
    fig.set_size_inches(7, 15)
    ax_list = ax_list.flat[:]
    for n, ax in enumerate(ax_list):
        in_force_comp = in_force_array[:, n]
        out_force_comp = out_force_array[:, n]

        ax.scatter(in_force_comp, out_force_comp)
        # get the appropriate limits for the plot
        for_limits = np.concatenate((in_force_comp, out_force_comp))
        flim = (for_limits.min() - 1, for_limits.max() + 1)
        ax.set_xlim(flim)
        ax.set_ylim(flim)
        # add line of
        ax.plot(flim, flim, c="k")
        # set labels
        ax.set_ylabel("Force by GAP [eV Å$^{-1}$]")
        ax.set_xlabel("Force by DFT [eV Å$^{-1}$]")
        # set title
        ax.set_title("Forces on " + symbol + ": " + components[n] + " components")
        # add text about RMSE
        _rms = rms_dict(in_force_comp, out_force_comp)
        rmse_text = (
            "RMSE:\n" + str(np.round(_rms["rmse"], 3)) + " eV Å$^{-1}$"
        )  # + ' +- ' + str(np.round(_rms['std'], 3)) + ' eV Å$^{-1}$'
        ax.text(
            0.9,
            0.1,
            rmse_text,
            transform=ax.transAxes,
            fontsize="large",
            horizontalalignment="right",
            verticalalignment="bottom",
        )

    plt.savefig(plotFolder + "Forces_SOAP_" + symbol + ".png", format="png", dpi=300)
    plt.show()


# fig, ax_list = plt.subplots(nrows=4, ncols=1, gridspec_kw={'hspace': 0.6})
# fig.set_size_inches(7, 20)
# ax_list = ax_list.flat[:]
plotFolder = sys.argv[1]
# plotFolder = "Iterative_Training_Figures_round1/"

# trainingData = "Training_Data" + "/" + "Iterative_Training_round1.xyz"
trainedData = "Trained_Data" + "/" + sys.argv[2]
trainingData = trainedData
energy_plot(trainingData, trainedData, "Energy of training data", plotFolder)
force_plot(trainingData, trainedData, "Si", "Forces on training data - Si", plotFolder)
force_plot(trainingData, trainedData, "H", "Forces on training data - H", plotFolder)
stress_plot(trainingData, trainedData, "Stress", plotFolder)
# plt.savefig('EnergiesAndForces_SOAP.png',format = 'png', dpi = 300)
# plt.show()
