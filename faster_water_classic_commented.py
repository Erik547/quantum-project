import time 
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.mappers import ParityMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_algorithms.optimizers import SLSQP, COBYLA, SciPyOptimizer
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_aer.primitives import Estimator
import matplotlib.pyplot as plt
import numpy as np
import psi4

from matplotlib import cbook, cm
from matplotlib.colors import LightSource

def classic_h_oh(memory='2 GB'):
    # Set memory and output file
    psi4.set_memory(memory)
    psi4.core.set_output_file('h2o_energy_plot.dat', False)
    tim = time.time()

    # Define a range of O-H bond distances to calculate the energy for
    # Note: Typically, the bond angle in water is 104.5 degrees, and we keep that constant here.
    # You might want to vary the bond angle as well for a more thorough analysis.

    X = np.linspace(0.98, 0.99, 20)
    fr = 98
    to = 105
    THETA = np.linspace(fr*np.pi/180, to*np.pi/180, 20)
    exact_energies = []
    min_dist = 100
    min_angle = 1000
    min_energy = 1000
    
    classic_energies = np.empty(shape=(len(X),len(THETA)))

    # Minimum distance: 1.331578947368421A Minimum angle: 3.087736779528254 (176.9 degrees) Energy:-15.759277375375431
    # Minimum distance: 1.3333333333333335A Minimum angle: 3.087736779528254 (176.9 degrees) Energy:-15.75928175414282
    # Minimum distance: 1.334A Minimum angle: 3.1380307571571517 (179.8 degrees) Energy:-15.759372006730736


    # Minimum distance: 0.9505263157894737A Minimum angle: 1.9198621771937625 (110.0 degrees) Energy:-75.98528111031735
    # Minimum distance: 0.9496666666666667A Minimum angle: 1.9428270357726352 (111.32 degrees) Energy:-75.98534244281039

    for (i, dist) in enumerate(X):
        for (j, angle) in enumerate(THETA):
            # Define the molecule at the current distance and bond angle
            h2o = psi4.geometry(f"""
                0 1
                O
                H 1 {dist}
                H 1 {dist} 2 {np.round(angle * 180/np.pi, 2)}
                """)

            # Set the basis and options
            psi4.set_options({'basis': 'sto-6g', 'reference': 'rhf'})

            # Calculate the energy
            energy = psi4.energy('scf')
            classic_energies[i][j] = energy
            exact_energies.append(energy)
            print(f"At distance {dist}A, Angle: {np.round(angle * 180/np.pi, 2)}, the energy is {energy}")
            if min_energy > energy:
                min_energy = energy
                min_angle = angle
                min_dist = dist

    # Set up plot
    print(f"Took {time.time() - tim}s classically")
    print(f"Minimum distance: {min_dist}A Minimum angle: {min_angle} ({np.round(min_angle * 180/np.pi, 2)} degrees) Energy:{min_energy}")
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

    ls = LightSource(270, 45)
    X, THETA = np.meshgrid(X, THETA)
    # To use a custom hillshading mode, override the built-in shading and pass
    # in the rgb colors of the shaded surface calculated from "shade".
    rgb = ls.shade(classic_energies.T, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
    surf = ax.plot_surface(X, THETA, classic_energies.T, rstride=1, cstride=1, facecolors=rgb,
                        linewidth=0, antialiased=False, shade=False)

    plt.xlabel('Interatomic distance (A)')
    plt.ylabel('Angle (rad)')
    plt.title('Energy (Hartree)')

    plt.show()
    # plt.plot(X, exact_energies, label='Exact Energy', color='blue')
    # plt.xlabel('Atomic distance (Angstrom)')
    # plt.ylabel('Energy (Hartree)')
    # plt.title("$Be$ ground state energy simulation")
    # plt.legend()
    # plt.show()

classic_h_oh()

# fr = 175
# to = 185

# X = np.linspace(1.20, 1.40, 50) # distance
# xlen = len(X)
# THETA = np.linspace(fr*np.pi/180, to*np.pi/180, 50) # angle in degrees
# theta_len = len(THETA)

# classic_energies = np.empty(shape=(xlen,theta_len))

# file1 = open('data.txt', 'r')

# for i in range(50):
#     for j in range(50):
#         classic_energies[i][j] = file1.readline()

# # f = open("data.txt", "r")
# # for (i, x) in enumerate(f):
# #     for (j, v) in enumerate(x):
# #         classic_energies[i][j] = v


# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

# ls = LightSource(270, 45)
# X, THETA = np.meshgrid(X, THETA)
# # To use a custom hillshading mode, override the built-in shading and pass
# # in the rgb colors of the shaded surface calculated from "shade".
# rgb = ls.shade(classic_energies.T, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
# surf = ax.plot_surface(X, THETA, classic_energies.T, rstride=1, cstride=1, facecolors=rgb,
#                     linewidth=0, antialiased=False, shade=False)

# ax.set_xlabel('Interatomic distance (A)', labelpad = 10)
# ax.set_ylabel('Angle (rad)', labelpad = 10)
# plt.title('Energy (Hartree)')
# ax.set_zlabel('Energy', labelpad = 20)
# # ax.set_zticks(ticks=20)
# ax.tick_params(axis='z', pad=10)

# plt.show()

# plt.figure(figsize=(8, 6))

# # Use 'extent' to specify the bounds of the x and y axes.
# im = plt.imshow(classic_energies.T, cmap=cm.gist_earth, aspect='auto',
#                 origin='lower', extent=[X.min(), X.max(), THETA.min(), THETA.max()])
# plt.colorbar(im)  # Add a colorbar to a plot
# plt.xlabel('Interatomic distance (A)')
# plt.ylabel('Angle (rad)')
# plt.title('Energy (Hartree)')
# plt.show()

    

# def get_qubit_op_H2O(angle, d = 0.957):
#     # Define the molecular structure of H2O
#     # Real O-H bond length in angstroms
    

#     # Define the coordinates
#     molecule = MoleculeInfo(
#         symbols=["H", "O", "H"],
#         coords=([d, 0, 0.0],  # Hydrogen
#                 [0,0,0],  # Oxygen at origin
#                 [d*np.cos(angle), d*np.sin(angle), 0.0]),
#         multiplicity=1,  # 2 * spin + 1; singlet state
#         charge=0,
#     ) #0.96
#     driver = PySCFDriver.from_molecule(molecule,basis='sto-6g')

#     # Get properties
#     properties = driver.run()

#     # Now you can get the reduced electronic structure problem
#     problem = FreezeCoreTransformer(
#         freeze_core=True, remove_orbitals=[-3, -2]
#     ).transform(properties)

#     num_particles = problem.num_particles
#     num_spatial_orbitals = problem.num_spatial_orbitals

#     mapper = BravyiKitaevMapper()
#     qubit_op = mapper.map(problem.second_q_ops()[0])
#     return qubit_op, num_particles, num_spatial_orbitals, problem, mapper



# angles = np.linspace(np.pi/2, 2*np.pi/3, 21)
# distances = np.linspace(0.75, 1.15, 41)
# # vqe_energies = []
# optimizer = COBYLA(maxiter=3)
# noiseless_estimator = Estimator(approximation=True)

# # atomic_distances, exact_energies = classic_h_oh()

# # for (i, angle) in enumerate(angles):
# #     tim = time.time()
# #     (qubit_op, num_particles, num_spatial_orbitals, problem, mapper) = get_qubit_op_H2O(angle, 0.9)

# #     init_state = HartreeFock(num_spatial_orbitals, num_particles, mapper)
# #     var_form = UCCSD(
# #         num_spatial_orbitals, num_particles, mapper, initial_state=init_state
# #     )
    
# #     vqe = VQE(
# #         noiseless_estimator,
# #         var_form,
# #         optimizer,
# #         initial_point=[0] * var_form.num_parameters,
# #     )
# #     tim2 = time.time()
# #     print("test", tim2-tim)
# #     vqe_calc = vqe.compute_minimum_eigenvalue(qubit_op)
# #     vqe_result = problem.interpret(vqe_calc).total_energies[0].real
# #     print("test2", time.time() - tim2)
# #     vqe_energies.append(vqe_result)
# #     print(f"Iteration: {i} Angle: {np.round(angle * 180/np.pi, 1)} VQE Result: {vqe_result:.5f}")

# # print("All energies have been calculated")


# # plt.plot(angles, vqe_energies,label = "VQE Energy")
# # #plt.plot(atomic_distances, exact_energies, label='Exact Energy', color='blue')
# # plt.xlabel('Atomic distance (Angstrom)')
# # plt.ylabel('Energy (Hartree)')
# # plt.title("$H_2O$ ground state energy simulation")
# # plt.legend()
# # plt.show()




# ax = plt.figure().add_subplot(projection='3d')

# # Make data.
# X = np.linspace(0.75, 1.15, 5) # distance
# xlen = len(X)
# THETA = np.linspace(np.pi/2, 2*np.pi/3, 6) # angle in degrees
# theta_len = len(THETA)

# vqe_energies = np.empty(shape=(xlen,theta_len))

# # Result

# for (i, x) in enumerate(X):
#     for (j, theta) in enumerate(THETA):
#         tim = time.time()
#         # print(x, theta)
#         (qubit_op, num_particles, num_spatial_orbitals, problem, mapper) = get_qubit_op_H2O(theta, x)

#         init_state = HartreeFock(num_spatial_orbitals, num_particles, mapper)
#         ansatz = UCCSD(
#             num_spatial_orbitals, num_particles, mapper, initial_state=init_state
#         )
        
#         vqe = VQE(
#             noiseless_estimator,
#             ansatz,
#             optimizer,
#             initial_point=[0] * ansatz.num_parameters,
#         )
#         vqe_calc = vqe.compute_minimum_eigenvalue(qubit_op)
#         vqe_result = problem.interpret(vqe_calc).total_energies[0].real
#         vqe_energies[i][j] = vqe_result
#         print(f"Interatomic distance: {x}A, Angle: {np.round(theta * 180/np.pi, 2)} VQE Energy: {vqe_result:.5f}, Time:{(time.time() - tim):.2f}s")

# X, THETA = np.meshgrid(X, THETA)

# # Set up plot
# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

# ls = LightSource(270, 45)
# # To use a custom hillshading mode, override the built-in shading and pass
# # in the rgb colors of the shaded surface calculated from "shade".
# rgb = ls.shade(vqe_energies, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
# surf = ax.plot_surface(X, THETA, vqe_energies.T, rstride=1, cstride=1, facecolors=rgb,
#                        linewidth=0, antialiased=False, shade=False)

# plt.show()
