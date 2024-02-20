The codes under this directory are used to generate the visualization files (.vtk) of the magnetic models and predict the near-Earth IMF parameters and geomagnetic index.

'SuModel2vtkSphere_Ou_Teng.pro' uses the IDL .sav file from the flux rope insertion method and convert it to a Stonyhurst spherical coordinates and save the data into the .vtk file, which could be visualized using Paraview.

'SuModel2vtkSphere_BottomBz.pro' convert the bottom magnetogram of the model into spherical coordinates and save to a .vtk file.

'fits2vtk_aia.pro' uses the AIA background image (.fits file) to generate .vtk file for visualization together with the magnetic model.

'PFSS_synoptic.py' uses the HMI synoptic magnetogram (.fits file) to extrapolate the global PFSS field. The magnetic field and decay index are contained in .dat files, which can be later transferred into .vtk files via 'PFSS_B_dat_2_vtk.pro' and 'PFSS_decay_index_dat_2_vtk.pro', respectively.

'SolarOrbiter_FluxRopeOptimizer.py' fits the in-situ observation by Solar Orbiter with the flux rope in the coronal magnetic field model, and then predict the near-Earth IMF parameters and geomagnetic index SYM-H (using Burton's Equations).

'SolarOrbiter_FluxRopeOptimizer_pseudo_prediction.py' are similar to 'SolarOrbiter_FluxRopeOptimizer.py' but uses the observed interplanetary IMF and solar wind parameters by WIND spacecraft to generate the SYM-H index, as shown in Figure 11.
