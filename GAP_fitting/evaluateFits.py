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


def evaluate_energies(in_file, out_file):
    """ Plots the distribution of energy per atom on the output vs the input"""
    # read files
    in_atoms = ase.io.read(in_file, ":")
    ener_in = [at.info["dft_energy"] / len(at) for at in in_atoms]
    ener_in_array = np.array(ener_in)
    # ener_in_array = np.delete(ener_in_array, ener_in_array.argmax())

    out_atoms = ase.io.read(out_file, ":")
    ener_out = [at.info["energy"] / len(at) for at in out_atoms]
    ener_out_array = np.array(ener_out)
    # ener_out_array = np.delete(ener_out_array, ener_out_array.argmax())

    n = 1
    e_rmse = []
    e_frames = []
    for energy_in, energy_out in zip(ener_in_array, ener_out_array):
        _rms = rms_dict(energy_in, energy_out)
        rmse_text = "RMSE: " + str(np.round(_rms["rmse"], 3)) + " eV/atom"
        # print("Frame ", n, " energy ", rmse_text)
        frame_text = "Frame " + str(n) + " energy " + rmse_text
        n += 1

        e_rmse.append(np.round(_rms["rmse"], 3))
        e_frames.append(frame_text)

    e_rmse_array = np.array(e_rmse)
    e_frames_array = np.array(e_frames)
    sorted_indices = np.argsort(e_rmse_array)[::-1]
    frames = e_frames_array[sorted_indices]

    for i in range(15):
        print(frames[i])


def evaluate_stresses(in_file, out_file):
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

    stress_rmse = []
    stress_frames = []
    for n in range(len(in_atoms)):
        stress_comp_in = stress_in_array[n, :]
        stress_comp_out = stress_out_array[n, :]
        _rms = rms_dict(stress_comp_in, stress_comp_out)
        rmse_text = "RMSE: " + str(np.round(_rms["rmse"], 3)) + " eV"
        frame_text = "Frame " + str(n + 1) + " stress " + rmse_text

        stress_rmse.append(np.round(_rms["rmse"], 3))
        stress_frames.append(frame_text)

    stress_rmse_array = np.array(stress_rmse)
    stress_frames_array = np.array(stress_frames)
    sorted_indices = np.argsort(stress_rmse_array)[::-1]
    frames = stress_frames_array[sorted_indices]

    for i in range(15):
        print(frames[i])

    # components = [
    #     r"$\sigma_{xx}$",
    #     r"$\sigma_{xy}$",
    #     r"$\sigma_{xz}$",
    #     r"$\sigma_{yx}$",
    #     r"$\sigma_{yy}$",
    #     r"$\sigma_{yz}$",
    #     r"$\sigma_{zx}$",
    #     r"$\sigma_{zy}$",
    #     r"$\sigma_{zz}$",
    # ]


def evaluate_forces(in_file, out_file, symbol="HSi"):
    """Plots the distribution of force components per atom on the output vs the input
    only plots for the given atom type(s)"""

    in_atoms_frames = [i for i in ase.io.iread(in_file, ":")]
    out_atoms_frames = [i for i in ase.io.iread(out_file, ":")]

    # extract data for only one species
    n = 1
    force_rmse = []
    force_frames = []
    for in_atoms, out_atoms in zip(in_atoms_frames, out_atoms_frames):
        in_force, out_force = [], []
        sym_all = in_atoms.get_chemical_symbols()
        for j, sym in enumerate(sym_all):
            if sym in symbol:
                # in_force.append(at_in.get_forces()[j])
                in_force.append(in_atoms.arrays["dft_force"][j])
                # out_force.append(at_out.get_forces()[j]) \
                out_force.append(
                    out_atoms.arrays["force"][j]
                )  # because QUIP and ASE use different names

        # inforce = in_force[n]
        # print("Length of in_force: ", len(in_force))
        # outforce = out_force[n]
        _rms = rms_dict(in_force, out_force)
        rmse_text = "RMSE: " + str(np.round(_rms["rmse"], 3)) + " eV/Å"
        frame_text = "Frame " + str(n) + " " + symbol + " force " + rmse_text
        n += 1

        force_rmse.append(np.round(_rms["rmse"], 3))
        force_frames.append(frame_text)

    force_rmse_array = np.array(force_rmse)
    force_frames_array = np.array(force_frames)
    sorted_indices = np.argsort(force_rmse_array)[::-1]
    frames = force_frames_array[sorted_indices]

    for i in range(15):
        print(frames[i])

    # in_force_array = np.asarray(in_force)
    # out_force_array = np.asarray(out_force)
    # components = [r"$\hat{x}$", r"$\hat{y}$", r"$\hat{z}$"]
    # for n, ax in enumerate(ax_list):
    #     in_force_comp = in_force_array[:, n]
    #     out_force_comp = out_force_array[:, n]


# trainingData = "Training_Data" + "/" + "Iterative_Training_round1.xyz"
trainedData = "Trained_Data" + "/" + sys.argv[1]
trainingData = trainedData
evaluate_energies(trainingData, trainedData)
evaluate_forces(trainingData, trainedData, "Si")
evaluate_forces(trainingData, trainedData, "H")
evaluate_stresses(trainingData, trainedData)
