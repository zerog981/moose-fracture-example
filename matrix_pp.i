# MATRIX - Benchmark
# C.Scherounigg
# Based on https://mooseframework.inl.gov/modules/porous_flow/multiapp_fracture_flow_introduction.html

conductive_heat_transfer = true
mass_transfer = true
enable_fractures = true
gravitational_acceleration = -9.81 # m/s²

water_fluid_properties = tabulated_water # tabulated_water, eos_water
gas_fluid_properties = tabulated_gas # tabulated_gas, eos_gas

[Mesh]
  uniform_refine = 0
  [matrix_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20
    xmin = -11
    xmax = 210
    ny = 20
    ymin = -10
    ymax = 160
    nz = 20
    zmin = -10
    zmax = 210
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [matrix_pressure_water]
    initial_condition = 28e6
  []
  [matrix_temperature]
    initial_condition = 673.15
    #scaling = 1E-6 # Fluid enthalpy is roughly 1E6
  []
  [matrix_pressure_gas]
    initial_condition = 28.1e6
  []
[]

[Kernels]
  [matrix_mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = matrix_pressure_water
  []
  [matrix_flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = matrix_pressure_water
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [matrix_mass_gas_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = matrix_pressure_gas
  []
  [matrix_flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = matrix_pressure_gas
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = matrix_temperature
  []
  [advection]
    type = PorousFlowHeatAdvection
    variable = matrix_temperature
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [conduction]
    type = PorousFlowHeatConduction
    variable = matrix_temperature
  []
[]

[DiracKernels]
  [heat_from_fracture]
    # Add heat received from fracture to matrix (Transfer Step 11a)
    type = ReporterPointSource
    variable = matrix_temperature
    value_name = heat_transfer_rate/transferred_heat_transfer_rate
    x_coord_name = heat_transfer_rate/x
    y_coord_name = heat_transfer_rate/y
    z_coord_name = heat_transfer_rate/z
    enable = ${conductive_heat_transfer}
  []
  [mass_water_from_fracture]
    # Add mass of water received from fracture to matrix (Transfer Step 11b)
    type = ReporterPointSource
    variable = matrix_pressure_water
    value_name = mass_flux_water/transferred_mass_flux_water
    x_coord_name = mass_flux_water/x
    y_coord_name = mass_flux_water/y
    z_coord_name = mass_flux_water/z
    enable = ${mass_transfer}
  []
  [mass_gas_from_fracture]
    # Add saturation of gas received from fracture to matrix (Transfer Step 11c)
    type = ReporterPointSource
    variable = matrix_pressure_gas
    value_name = mass_flux_gas/transferred_mass_flux_gas
    x_coord_name = mass_flux_gas/x
    y_coord_name = mass_flux_gas/y
    z_coord_name = mass_flux_gas/z
    enable = ${mass_transfer}
  []
[]

[AuxVariables]
  [matrix_massfraction_ph0_sp0] # Phase 0, Species 0
    initial_condition = 1
  []
  [matrix_massfraction_ph1_sp0] # Phase 1, Species 0
    initial_condition = 0
  []
  [matrix_massfraction_ph0_sp1] # Phase 0, Species 1
  []
  [matrix_massfraction_ph1_sp1] # Phase 1, Species 1
  []
  [matrix_saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_saturation_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_density_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_density_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_enthalpy_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_thermal_conductivity]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_permeability]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_normal_x]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [fracture_normal_y]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  []
  [fracture_normal_z]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [element_normal_length]
    order = CONSTANT
    family = MONOMIAL
  []
  #[relative_permeability_material_water_aux]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[]
  #[rho_k_krel_over_mu]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[]
[]

[AuxKernels]
  [matrix_saturation_gas]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 1
    variable = matrix_saturation_gas
  []
  [matrix_saturation_water]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = matrix_saturation_water
  []
  [matrix_enthalpy_gas]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 1
    variable = matrix_enthalpy_gas
  []
  [matrix_enthalpy_water]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 0
    variable = matrix_enthalpy_water
  []
  [matrix_density_water]
    type = PorousFlowPropertyAux
    property = density
    phase = 0
    variable = matrix_density_water
  []
  [matrix_density_gas]
    type = PorousFlowPropertyAux
    property = density
    phase = 1
    variable = matrix_density_gas
  []
  [matrix_viscosity_water]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 0
    variable = matrix_viscosity_water
  []
  [matrix_viscosity_gas]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 1
    variable = matrix_viscosity_gas
  []
  [element_normal_length_auxkernel]
    # Calculate normal length in direction of fracture normal (Transfer Step 3) 
    type = PorousFlowElementLength
    variable = element_normal_length
    direction = 'fracture_normal_x fracture_normal_y fracture_normal_z'
  []
  [normal_thermal_conductivity_auxkernel]
    # Thermal conductivity in normal direction to fracture. Assumed constant instead of deriving it from a thermal conductivity tensor. (Transfer Step 4a)
    type = ConstantAux
    variable = normal_thermal_conductivity
    value = 5
    enable = ${conductive_heat_transfer}
  []
  [normal_permeability]
    # Permeability in normal direction to fracture. Assumed constant instead of deriving it from a permeability tensor. (Transfer Step 4b)
    type = ConstantAux
    variable = normal_permeability
    value = 1e-15
    enable = ${mass_transfer}
  []
  #[relative_permeability_water_material_auxkernel]
  #  # Use relative permeability material property as AuxVariable
  #  type = MaterialRealAux
  #  property = PorousFlow_relative_permeability_qp0
  #  variable = relative_permeability_material_water_aux
  #[]
  #[rho_k_krel_over_mu_auxkernel]
  #  type = ParsedAux
  #  variable = rho_k_krel_over_mu
  #  coupled_variables = 'matrix_density_water matrix_viscosity_water relative_permeability_material_water_aux normal_permeability'
  #  expression = 'matrix_density_water*normal_permeability*relative_permeability_material_water_aux/matrix_viscosity_water'
  #  enable = ${mass_transfer}
  #[]
[]

[Materials]
  [matrix_temperature]
    type = PorousFlowTemperature
    temperature = matrix_temperature
  []
  [pore_pressures]
    type = PorousFlow2PhasePP
    phase0_porepressure = matrix_pressure_water
    phase1_porepressure = matrix_pressure_gas
    capillary_pressure = capillary_pressure
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'matrix_massfraction_ph0_sp0 matrix_massfraction_ph1_sp0'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = ${water_fluid_properties} # Fluid properties
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = ${gas_fluid_properties} # Fluid properties
    phase = 1
  []
  [relative_permeability_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
    s_res = 0.1
    sum_s_res = 0.1
  []
  [relative_permeability_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
  []
  [porosity]
    type = PorousFlowPorosityConst # To supply initial value
    porosity = 0.01
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1E-10
    fluid_bulk_modulus = 2E9
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-15 0 0   0 1e-15 0   0 0 1e-15'
  []
  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    fluid_coefficient = 5E-6
    drained_coefficient = 2E-4
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1 0 0  0 1 0  0 0 1'
    wet_thermal_conductivity = '3 0 0  0 3 0  0 0 3'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
  []
[]

[FluidProperties]
  [eos_water]
    type = Water97FluidProperties
  []
  [eos_gas]
    type = CO2FluidProperties
  []
  [tabulated_water]
    type = TabulatedBicubicFluidProperties
    fp = eos_water
    fluid_property_file = fluid_properties_water.csv
    interpolated_properties = 'density enthalpy internal_energy viscosity'
    temperature_min = 300 # K
    temperature_max = 1000 # K
    pressure_min = 1e6 # Pa
    pressure_max = 50e6 # Pa
    num_T = 100
    num_p = 100
    error_on_out_of_bounds = False
  []
  [tabulated_gas]
    type = TabulatedBicubicFluidProperties
    fp = eos_gas
    fluid_property_file = fluid_properties_gas.csv
    interpolated_properties = 'density enthalpy internal_energy viscosity'
    temperature_min = 300 # K
    temperature_max = 1000 # K
    pressure_min = 1e6 # Pa
    pressure_max = 50e6 # Pa
    num_T = 100
    num_p = 100
    error_on_out_of_bounds = False
  []
  [simple_water]
    type = SimpleFluidProperties
  []
  [simple_gas]
    type = SimpleFluidProperties
    cp = 816.9
    cv = 627.9
    density0 = 1.997
    molar_mass = 0.018015
    thermal_conductivity = 0.0147
    viscosity = 13.72e-6
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'matrix_pressure_water matrix_temperature matrix_pressure_gas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [capillary_pressure]
    type = PorousFlowCapillaryPressureVG
    alpha = 1e-6
    m = 0.6
    s_scale = 0.9
  []
[]

[BCs]
  [side_pressure]
    type = DirichletBC
    variable = matrix_pressure_water
    value = 28e6
    boundary = 'front' # Corresponds to top of reservoir
  []
  [side_temperature]
    type = DirichletBC
    variable = matrix_temperature
    value = 673.15
    boundary = 'back' # Corresponds to top and bottom of reservoir
  []
[]

[Postprocessors]

[]

[VectorPostprocessors]
  [heat_transfer_rate]
    type = ConstantVectorPostprocessor
    vector_names = 'transferred_heat_transfer_rate x y z'
    value = '0; 0; 0; 0'
    outputs = none
    enable = ${conductive_heat_transfer}
  []
  [mass_flux_water]
    type = ConstantVectorPostprocessor
    vector_names = 'transferred_mass_flux_water x y z'
    value = '0; 0; 0; 0'
    outputs = none
    enable = ${mass_transfer}
  []
  [mass_flux_gas]
    type = ConstantVectorPostprocessor
    vector_names = 'transferred_mass_flux_gas x y z'
    value = '0; 0; 0; 0'
    outputs = none
    enable = ${mass_transfer}
  []
[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 9.4608e8 # 30 years
  nl_max_its = 20 # Maximum number of non-linear iterations
  l_max_its = 1000
  nl_abs_tol = 1e-6
  automatic_scaling = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 4
    growth_factor = 1.1
    cutback_factor = 0.5
  []
[]

[Outputs]
  print_linear_residuals = false
  checkpoint = false # Add checkpoint files for recovering simulation
  exodus = true
  [csv]
    type = CSV
  []
[]

[MultiApps]
  [fracture_app]
    type = TransientMultiApp
    input_files = fracture_pp.i
    cli_args = 'Outputs/ex/sync_only=false'
    execute_on = TIMESTEP_BEGIN
    sub_cycling = true
    enable = ${enable_fractures}
  []
[]

[Transfers]
  [element_normal_length_to_fracture]
    # Transfer length along fracture normal to fracture network (Transfer Step 5a)
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = fracture_app
    source_variable = element_normal_length
    variable = enclosing_element_normal_length
  []
  [element_normal_thermal_conductivity_to_fracture]
    # Transfer thermal conductivity along fracture normal to fracture network (Transfer Step 5b)
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = fracture_app
    source_variable = normal_thermal_conductivity
    variable = enclosing_element_normal_thermal_conductivity
    enable = ${conductive_heat_transfer}
  []
  [element_normal_permeability_to_fracture]
    # Transfer permeability along fracture normal to fracture network (Transfer Step 5c)
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = fracture_app
    source_variable = normal_permeability
    variable = enclosing_element_normal_permeability
    enable = ${mass_transfer}
  []
  #[rho_k_krel_over_mu_to_fracture]
  #  # Transfer expression for mass transfer coefficient along fracture normal to fracture network (Transfer Step 5c)
  #  type = MultiAppGeneralFieldNearestLocationTransfer
  #  to_multi_app = fracture_app
  #  source_variable = rho_k_krel_over_mu
  #  variable = enclosing_element_normal_rho_k_krel_over_mu
  #  enable = ${mass_transfer}
  #[]
  [temperature_to_fracture]
    # Transfer matrix temperature to fracture network. (Transfer Step 7a)
    type = MultiAppGeometricInterpolationTransfer
    to_multi_app = fracture_app
    source_variable = matrix_temperature
    variable = matrix_temperature_transfer
    enable = ${conductive_heat_transfer}
  []
  [pressure_water_to_fracture]
    # Transfer water pressure to fracture network (Transfer Step 7b)
    type = MultiAppGeometricInterpolationTransfer
    to_multi_app = fracture_app
    source_variable = matrix_pressure_water
    variable = matrix_pressure_water_transfer
    enable = ${mass_transfer}
  []
  [pressure_gas_to_fracture]
    # Transfer water pressure to fracture network (Transfer Step 7b)
    type = MultiAppGeometricInterpolationTransfer
    to_multi_app = fracture_app
    source_variable = matrix_pressure_gas
    variable = matrix_pressure_gas_transfer
    enable = ${mass_transfer}
  []
  [normal_x_from_fracture]
    # Transfer normal vector of fracture from fracture network (Transfer Step 2)
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = fracture_app
    source_variable = normal_dir_x
    variable = fracture_normal_x
  []
  [normal_y_from_fracture]
    # Transfer normal vector of fracture from fracture network (Transfer Step 2)
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = fracture_app
    source_variable = normal_dir_y
    variable = fracture_normal_y
  []
  [normal_z_from_fracture]
    # Transfer normal vector of fracture from fracture network (Transfer Step 2)
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = fracture_app
    source_variable = normal_dir_z
    variable = fracture_normal_z
  []
  [heat_from_fracture]
    # Receive heat from fracture network (Transfer Step 10a)
    type = MultiAppReporterTransfer
    from_multi_app = fracture_app
    from_reporters = 'heat_transfer_rate/heat_transfer_rate heat_transfer_rate/x heat_transfer_rate/y heat_transfer_rate/z'
    to_reporters = 'heat_transfer_rate/transferred_heat_transfer_rate heat_transfer_rate/x heat_transfer_rate/y heat_transfer_rate/z'
    enable = ${conductive_heat_transfer}
  []
  [mass_water_from_fracture]
    # Receive mass of water from fracture network (Transfer Step 10b)
    type = MultiAppReporterTransfer
    from_multi_app = fracture_app
    from_reporters = 'mass_flux_water/mass_flux_water mass_flux_water/x mass_flux_water/y mass_flux_water/z'
    to_reporters = 'mass_flux_water/transferred_mass_flux_water mass_flux_water/x mass_flux_water/y mass_flux_water/z'
    enable = ${mass_transfer}
  []
  [mass_gas_from_fracture]
    # Receive mass of water from fracture network (Transfer Step 10c)
    type = MultiAppReporterTransfer
    from_multi_app = fracture_app
    from_reporters = 'mass_flux_gas/mass_flux_gas mass_flux_gas/x mass_flux_gas/y mass_flux_gas/z'
    to_reporters = 'mass_flux_gas/transferred_mass_flux_gas mass_flux_gas/x mass_flux_gas/y mass_flux_gas/z'
    enable = ${mass_transfer}
  []
[]
