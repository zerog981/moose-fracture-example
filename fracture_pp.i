# FRACTURE
# C.Scherounigg
# Based on https://mooseframework.inl.gov/modules/porous_flow/multiapp_fracture_flow_introduction.html

conductive_heat_transfer = true
mass_transfer = true
enable_production = true
gravitational_acceleration = -9.81 # m/s²
injection_rate = 15 # kg/s
injection_temperature = 673.15 # K
production_pressure = 23e6 # Pa
water_fluid_properties = tabulated_water # tabulated_water, eos_water, simple_water
gas_fluid_properties = tabulated_gas # tabulated_gas, eos_gas, simple_gas
[Mesh]
  uniform_refine = 0
  [test_fracture_mesh]
    type = FileMeshGenerator
    file = 'Cluster_34.exo'
  []
  [injection_node]
    type = BoundingBoxNodeSetGenerator
    input = test_fracture_mesh
    bottom_left = '-1000 0 -1000'
    top_right = '1000 0.504 1000'
    new_boundary = injection_node
  []
[]
[GlobalParams]
  PorousFlowDictator = dictator
[]
[Functions]

[]
[Variables]
  [fracture_pressure_water]
    initial_condition = 28e6
  []
  [fracture_temperature]
    initial_condition = 673.15
  []
  [fracture_pressure_gas]
    initial_condition = 28.1e6
  []
[]
[Kernels]
  [fracture_mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = fracture_pressure_water
  []
  [fracture_flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = fracture_pressure_water
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [fracture_mass_gas_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = fracture_pressure_gas
  []
  [fracture_flux_gas]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = fracture_pressure_gas
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = fracture_temperature
  []
  [advection]
    type = PorousFlowHeatAdvection
    variable = fracture_temperature
    gravity = '0 0 ${gravitational_acceleration}'
  []
  [conduction]
    type = PorousFlowHeatConduction
    variable = fracture_temperature
  []
  [heat_to_matrix]
    # Transfer heat to matrix (Transfer Step 8a)
    type = PorousFlowHeatMassTransfer
    variable = fracture_temperature
    v = matrix_temperature_transfer
    transfer_coefficient = heat_transfer_coefficient
    save_in = heat_transfer_rate # Heat transfer rate in W.
    enable = ${conductive_heat_transfer}
  []
  [water_to_matrix]
    # Transfer mass of water to matrix (Transfer Step 8b)
    type = PorousFlowHeatMassTransfer
    variable = fracture_pressure_water
    v = matrix_pressure_water_transfer
    transfer_coefficient = mass_transfer_coefficient_water
    save_in = mass_flux_water # Mass flux in kg/(m²s).
    enable = ${mass_transfer}
  []
  [gas_to_matrix]
    # Transfer mass of gas to matrix (Transfer Step 8c)
    type = PorousFlowHeatMassTransfer
    variable = fracture_pressure_gas
    v = matrix_pressure_gas_transfer
    transfer_coefficient = mass_transfer_coefficient_gas
    save_in = mass_flux_gas # Mass flux in kg/(m²s).
    enable = ${mass_transfer}
  []
[]
[DiracKernels]
  [injection]
    type = PorousFlowSquarePulsePointSource
    mass_flux = ${injection_rate}
    point = '58.8124 0.50384 74.7838'
    variable = fracture_pressure_gas
  []
  [water_production]
    type = PorousFlowPeacemanBorehole
    SumQuantityUO = mass_water_produced
    variable = fracture_pressure_water
    function_of = pressure
    bottom_p_or_t = ${production_pressure}
    character = 1
    line_length = 1
    point_file = production.xyz
    unit_weight = '0 0 0'
    fluid_phase = 0
    use_mobility = true
    enable = ${enable_production}
  []
  [gas_production]
    type = PorousFlowPeacemanBorehole
    SumQuantityUO = mass_gas_produced
    variable = fracture_pressure_gas
    function_of = pressure
    bottom_p_or_t = ${production_pressure}
    character = 1
    line_length = 1
    point_file = production.xyz
    unit_weight = '0 0 0'
    fluid_phase = 1
    use_mobility = true
    enable = ${enable_production}
  []
  [heat_production_water]
    type = PorousFlowPeacemanBorehole
    SumQuantityUO = heat_water_produced
    variable = fracture_temperature
    function_of = pressure
    bottom_p_or_t = ${production_pressure}
    character = 1
    line_length = 1
    point_file = production.xyz
    unit_weight = '0 0 0'
    fluid_phase = 0
    use_mobility = true
    use_enthalpy = true
    enable = ${enable_production}
  []
  [heat_production_gas]
    type = PorousFlowPeacemanBorehole
    SumQuantityUO = heat_gas_produced
    variable = fracture_temperature
    function_of = pressure
    bottom_p_or_t = ${production_pressure}
    character = 1
    line_length = 1
    point_file = production.xyz
    unit_weight = '0 0 0'
    fluid_phase = 1
    use_mobility = true
    use_enthalpy = true
    enable = ${enable_production}
  []
[]
[AuxVariables]
  [fracture_massfraction_ph0_sp0]
    initial_condition = 1
  []
  [fracture_massfraction_ph1_sp0]
    initial_condition = 0
  []
  [fracture_massfraction_ph0_sp1]
  []
  [fracture_massfraction_ph1_sp1]
  []
  [fracture_saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_saturation_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_density_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_density_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [fracture_enthalpy_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [matrix_temperature_transfer]
    initial_condition = 673.15
    enable = ${conductive_heat_transfer}
  []
  [matrix_pressure_water_transfer]
    initial_condition = 28e6
    enable = ${mass_transfer}
  []
  [matrix_pressure_gas_transfer]
    initial_condition = 28.1e6
    enable = ${mass_transfer}
  []
  [heat_transfer_coefficient]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
    enable = ${conductive_heat_transfer}
  []
  [mass_transfer_coefficient_water]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
    enable = ${mass_transfer}
  []
  [mass_transfer_coefficient_gas]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
    enable = ${mass_transfer}
  []
  [enclosing_element_normal_length]
    order = CONSTANT
    family = MONOMIAL
  []
  [enclosing_element_normal_thermal_conductivity]
    order = CONSTANT
    family = MONOMIAL
  []
  [enclosing_element_normal_permeability]
    order = CONSTANT
    family = MONOMIAL
  []
#  [enclosing_element_normal_rho_k_krel_over_mu]
#    order = CONSTANT
#    family = MONOMIAL
#  []
  [heat_transfer_rate]
    initial_condition = 0.0
  []
  [mass_flux_water]
    initial_condition = 0.0
  []
  [mass_flux_gas]
    initial_condition = 0.0
  []
  [normal_dir_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_dir_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_dir_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [rho_water_times_gravity]
    order = CONSTANT
    family = MONOMIAL
  []
  [relative_permeability_water_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [relative_permeability_gas_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [normal_dir_x_auxkernel]
    # Calculate normal vector of fracture (Transfer Step 1)
    type = PorousFlowElementNormal
    variable = normal_dir_x
    component = x
  []
  [normal_dir_y_auxkernel]
    # Calculate normal vector of fracture (Transfer Step 1)
    type = PorousFlowElementNormal
    variable = normal_dir_y
    component = y
  []
  [normal_dir_z_auxkernel]
    # Calculate normal vector of fracture (Transfer Step 1)
    type = PorousFlowElementNormal
    variable = normal_dir_z
    component = z
  []
  [fracture_saturation_gas]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 1
    variable = fracture_saturation_gas
  []
  [fracture_saturation_water]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = fracture_saturation_water
  []
  [fracture_enthalpy_gas]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 1
    variable = fracture_enthalpy_gas
  []
  [fracture_enthalpy_water]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 0
    variable = fracture_enthalpy_water
  []
  [fracture_density_water]
    type = PorousFlowPropertyAux
    property = density
    phase = 0
    variable = fracture_density_water
  []
  [fracture_density_gas]
    type = PorousFlowPropertyAux
    property = density
    phase = 1
    variable = fracture_density_gas
  []
  [fracture_viscosity_water]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 0
    variable = fracture_viscosity_water
  []
  [fracture_viscosity_gas]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 1
    variable = fracture_viscosity_gas
  []
  [heat_transfer_coefficient_auxkernel]
    # Calculate heat transfer coefficient (Transfer Step 6a)
    type = ParsedAux
    variable = heat_transfer_coefficient
    coupled_variables = 'enclosing_element_normal_length enclosing_element_normal_thermal_conductivity'
    constant_names = h_s
    constant_expressions = 1e3
    expression = 'if(enclosing_element_normal_length = 0, 0, h_s * enclosing_element_normal_thermal_conductivity * 2 * enclosing_element_normal_length / (h_s * enclosing_element_normal_length * enclosing_element_normal_length * enclosing_element_normal_thermal_conductivity * 2 * enclosing_element_normal_length))'
    enable = ${conductive_heat_transfer}
  []
  [mass_transfer_coefficient_water_auxkernel]
    # Calculate mass transfer coefficient for water (Transfer Step 6b)
    type = ParsedAux
    variable = mass_transfer_coefficient_water
    coupled_variables = 'enclosing_element_normal_length enclosing_element_normal_permeability fracture_density_water fracture_viscosity_water relative_permeability_water_aux'
    expression = 'if(enclosing_element_normal_length = 0, 0, 2 * fracture_density_water * enclosing_element_normal_permeability * relative_permeability_water_aux / (fracture_viscosity_water * enclosing_element_normal_length))'
    #coupled_variables = 'enclosing_element_normal_rho_k_krel_over_mu enclosing_element_normal_length'
    #expression = 'if(enclosing_element_normal_length = 0, 0, 2 / enclosing_element_normal_length * enclosing_element_normal_rho_k_krel_over_mu)'
    enable = ${mass_transfer}
  []
  [mass_transfer_coefficient_gas_auxkernel]
    # Calculate mass transfer coefficient for gas (Transfer Step 6c)
    type = ParsedAux
    variable = mass_transfer_coefficient_gas
    coupled_variables = 'enclosing_element_normal_length enclosing_element_normal_permeability fracture_density_gas fracture_viscosity_gas relative_permeability_gas_aux'
    expression = 'if(enclosing_element_normal_length = 0, 0, 2 * fracture_density_gas * enclosing_element_normal_permeability * relative_permeability_gas_aux / (fracture_viscosity_gas * enclosing_element_normal_length))'
    enable = ${mass_transfer}
  []
  [relative_permeability_water_material_aux]
    # Use relative permeability material property as AuxVariable
    type = MaterialRealAux
    property = PorousFlow_relative_permeability_qp0
    variable = relative_permeability_water_aux
  []
  [relative_permeability_gas_material_aux]
    # Use relative permeability material property as AuxVariable
    type = MaterialRealAux
    property = PorousFlow_relative_permeability_qp1
    variable = relative_permeability_gas_aux
  []
[]
[Materials]
  [fracture_temperature]
    type = PorousFlowTemperature
    temperature = fracture_temperature
  []
  [pore_pressures]
    type = PorousFlow2PhasePP
    phase0_porepressure = fracture_pressure_water
    phase1_porepressure = fracture_pressure_gas
    capillary_pressure = capillary_pressure
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'fracture_massfraction_ph0_sp0 fracture_massfraction_ph1_sp0'
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
    porosity = 1e-1
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1E-10
    fluid_bulk_modulus = 2E9
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0   0 1e-12 0   0 0 1e-12'
  []
  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    fluid_coefficient = 5E-6
    drained_coefficient = 2E-4
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0.6e-4 0 0  0 0.6e-4 0  0 0 0.6e-4'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 0
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
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
  [mass_water_produced]
    type = PorousFlowSumQuantity
  []
  [mass_gas_produced]
    type = PorousFlowSumQuantity
  []
  [heat_water_produced]
    type = PorousFlowSumQuantity
  []
  [heat_gas_produced]
    type = PorousFlowSumQuantity
  []
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'fracture_pressure_water fracture_temperature fracture_pressure_gas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [capillary_pressure]
    #type = PorousFlowCapillaryPressureConst
    #pc = 0
    #type = PorousFlowCapillaryPressureBC
    #lambda = 1
    #pe = 50000
    type = PorousFlowCapillaryPressureVG
    alpha = 1e-6
    m = 0.6
    s_scale = 0.9
  []
[]
[BCs]
  [injection_temperature]
    type = DirichletBC
    variable = fracture_temperature
    value = ${injection_temperature}
    boundary = 'injection_node'
  []
[]
[ICs]
  [injection_temperature]
    type = ConstantIC
    variable = fracture_temperature
    value = ${injection_temperature}
    boundary = 'injection_node'
  []
[]
[Postprocessors]
  [injection_saturation_gas]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_saturation_gas
  []
  [injection_density_gas]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_density_gas
  []
  [injection_density_water]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_density_water
  []
  [injection_temperature]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_temperature
  []
  [injection_viscosity_water]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_viscosity_water
  []
  [injection_viscosity_gas]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_viscosity_gas
  []
  [injection_enthalpy_water]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_enthalpy_water
  []
  [injection_enthalpy_gas]
    type = PointValue
    point = '58.8124 0.50384 74.7838'
    variable = fracture_enthalpy_gas
  []
[]
[VectorPostprocessors]
  # Send heat to matrix (Transfer Step 9a)
  [heat_transfer_rate]
    type = NodalValueSampler
    outputs = none
    sort_by = id
    variable = heat_transfer_rate
    enable = ${conductive_heat_transfer}
  []
  # Send mass of water to matrix (Transfer Step 9b)
  [mass_flux_water]
    type = NodalValueSampler
    outputs = none
    sort_by = id
    variable = mass_flux_water
    enable = ${mass_transfer}
  []
  # Send mass of gas to matrix (Transfer Step 9b)
  [mass_flux_gas]
    type = NodalValueSampler
    outputs = none
    sort_by = id
    variable = mass_flux_gas
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
  nl_max_its = 15 # Maximum number of non-linear iterations
  l_max_its = 1000
  nl_abs_tol = 1e-6
  automatic_scaling = true
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.5
  []
[]
[Outputs]
  print_linear_residuals = false
  checkpoint = false # Create checkpoint files for recovering simulation
  [ex]
    type = Exodus
  []
  [csv]
    type = CSV
  []
[]
[Debug]
  show_material_props = false
[]
