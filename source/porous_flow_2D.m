%% Equation set
%
% Conservation of Momentum
%   d/dx( K * dP/dx ) = 0       v = - K * grad(P)
%
% Conservation of Energy
%   dT/dt - d/dx( (eb*Kw + (1-eb)*Ks)/(eb*rho*cpw + (1-eb)*rhos*cps) * dT/dx ) 
%           - eb*rho*cpw/(eb*rho*cpw + (1-eb)*rhos*cps)*v*dT/dx 
%           + eb/(eb*rho*cpw + (1-eb)*rhos*cps)*(-dHr)*[rxn (mol/m^3/s)] 
%
% Conservation of Mass
%   dCi/dt - d/dx ( m*Deff * dCi/dx ) = - m*v*dCi/dx + [rxn (mol/m^3/s)] 
%

% Mesh and Geometry files 
% https://www.mathworks.com/help/pde/geometry-and-mesh.html

classdef porous_flow_2D < handle 

    %% Public properties
    properties (Access = public)
        particle_diameter     % particle diameter in m
        particle_porosity     % Volume of micro-voids per bulk solids volume
        muA                   % Viscosity pre-exponential in kg/m/s
        muB                   % Viscosity temperature relation B in K
        muC                   % Viscosity temperature relation C in K 
        rhoA                  % Density temperature relation A in K^-2
        rhoB                  % Density temperature relation B in K^-1
        rhoC                  % Density temperature relation C in -
        rhoRef                % Reference Density in kg/m^3
        bulk_porosity         % system bulk porosity (vol voids / vol bulk solids)
        KozenyCarmannConst    % Kozeny-Carmann coefficient (dimensionless)
        diff_factor           % Diffusivity factor 
        alpha                 % Dispersion length factor in m
        refDiff               % Reference diffusivity in m^2/s
        refDiffTemp           % Reference diffusivity temperature in K
        cpA                   % Specific heat temperature parameter in K^-2
        cpB                   % Specific heat temperature parameter in K^-1
        cpC                   % Specific heat reference parameter -
        KwA                   % Thermal conductivity of water slope in K^-1
        KwB                   % Thermal conductivity of water intercept -
        char_len              % A characteristic length for the domain in m
        bulk_solids_dens      % Bulk solids density in the domain in kg/m^3
        bulk_solids_Cp        % Bulk solids heat capacity in J/kg/K
        bulk_solids_cond      % Bulk solids thermal conductivity in J/m/K/s

        mobile_spec_idx       % Matrix of 1s and 0s (where 1 indicates the species is mobile)
        rxn_stoich            % Stoichiometry matrix
        rxn_powers            % Matrix for reaction power terms
        rxn_rate_const        % Matrix for reaction rate constants (units will vary)
        rxn_act_energy        % Matrix for reaction rate activation energies (J/mol)
        rxn_enthalpy          % Matrix for reaction enthalpies (J/mol)

        model                 % PDE model object 
    end
    
    %% Public methods
    methods (Access = public)

        %% Constructor
        function obj = porous_flow_2D(N, Rxns, subdomains)
            % Validate the inputs
            arguments
                N (1,1) {mustBePositive} = 1
                Rxns (1,1) {mustBeNonnegative} = 0
                subdomains (1,1) {mustBePositive} = 1
            end
            % Set defaults 
            obj.bulk_porosity = ones(1,1,subdomains)*0.5;        % -
            obj.particle_diameter = ones(1,1,subdomains)*5e-3;   % m
            obj.muA = ones(1,1,subdomains)*0.02939/1000;         % kg/m/s
            obj.muB = ones(1,1,subdomains)*507.88;               % K
            obj.muC = ones(1,1,subdomains)*149.3;                % K
            obj.rhoA = ones(1,1,subdomains)*-2.9335E-6;          % K^-2
            obj.rhoB = ones(1,1,subdomains)*0.001529811;         % K^-1
            obj.rhoC = ones(1,1,subdomains)*0.787973;            % dimensionless
            obj.rhoRef = ones(1,1,subdomains)*1000;              % kg/m^3
            obj.KozenyCarmannConst = ones(1,1,subdomains)*5.55;  % - 
            obj.diff_factor = ones(1,1,subdomains)*0.5;          % -
            obj.alpha = ones(1,1,subdomains)*0.01;               % m
            obj.cpA = ones(1,1,subdomains)*0.01356461398877;     % K^-2
            obj.cpB = ones(1,1,subdomains)*-8.75598978135192;    % K^-1
            obj.cpC = ones(1,1,subdomains)*5591.8143667457;      % -
            obj.KwA = ones(1,1,subdomains)*0.001187275697628;    % K^-1
            obj.KwB = ones(1,1,subdomains)*0.247611225358853;    % - 
            obj.char_len = ones(1,1,subdomains)*5e-3;            % m
            obj.bulk_solids_dens = ones(1,1,subdomains)*2600;    % kg/m^3
            obj.bulk_solids_Cp = ones(1,1,subdomains)*1.315e3;   % J/kg/K
            obj.bulk_solids_cond = ones(1,1,subdomains)*1.04;    % J/m/K/s

            % These can be different for each chemical 
            obj.refDiff = ones(N,1,subdomains)*2.296E-5/100^2;   % m^2/s
            obj.refDiffTemp = ones(N,1,subdomains)*298.15;       % K

            % Setup space for reactions
            obj.mobile_spec_idx = ones(N,1,subdomains);
            obj.rxn_stoich = zeros(N,Rxns,subdomains);
            obj.rxn_powers = zeros(N,Rxns,subdomains);
            obj.rxn_rate_const = zeros(1,Rxns,subdomains);
            obj.rxn_act_energy = zeros(1,Rxns,subdomains);
            obj.rxn_enthalpy = zeros(1,Rxns,subdomains);

            obj.model = createpde(N+2);
        end

        %% Function to set model geometry from edges
        function [] = set_geometry_from_edges(obj, geo)
            geometryFromEdges(obj.model,geo);
        end

        
        %% Viscosity function of water
        %   @param temperature in K
        function vis = ViscosityWater(obj, temperature, sub)
            % Validate the inputs
            arguments
                obj
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            vis = obj.muA(1,1,sub)*exp(obj.muB(1,1,sub)./(temperature - obj.muC(1,1,sub)));
        end

        %% Density function of water
        %   @param pressure in Pa
        %   @param temperature in K
        function dens = DensityWater(obj, pressure, temperature, sub)
            % Validate the inputs
            arguments
                obj
                pressure {mustBeNumeric} = 101350  
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            a = 1.0135 + 4.9582e-7 * pressure;
            coeff = obj.rhoA(1,1,sub)*temperature.*temperature + obj.rhoB(1,1,sub)*temperature + obj.rhoC(1,1,sub);
            dens = a.*coeff*obj.rhoRef(1,1,sub);
        end

        %% Wall heat transfer coefficient for water (J/m^2/K/s)
        %
        %       Nu = 0.332 * Re_x^(1/2)*Pr^(1/3)
        %
        %               Nu = h*L/Kw
        %               Pr = cp*mu/Kw
        %               Re = rho*u_mag*L/mu
        %
        %   @param pressure in Pa
        %   @param temperature in K
        function h = WallHeatTransferWater(obj, pressure, temperature, ux, uy, uz, sub)
            % Validate the inputs
            arguments
                obj
                pressure {mustBeNumeric} = 101350  
                temperature {mustBeNumeric} = 298
                ux {mustBeNumeric} = 0
                uy {mustBeNumeric} = 0
                uz {mustBeNumeric} = 0
                sub {mustBePositive} = 1
            end
            rho = obj.DensityWater(pressure,temperature,sub);
            mu = obj.ViscosityWater(temperature,sub);
            Kw = obj.ThermalConductivityWater(temperature,sub);
            cp = obj.SpecificHeatWater(temperature,sub);
            umag = sqrt(ux.^2 + uy.^2 + uz.^2);
            Re = rho.*umag*obj.char_len(1,1,sub)./mu;
            Pr = cp.*mu./Kw;
            Nu = 3 + 0.332*Re.^(1/2).*Pr.^(1/3);
            h = Nu.*Kw/obj.char_len(1,1,sub);
        end
        
        %% Thermal conductivity function of water (in J/m/K/s)
        %   @param temperature in K
        function Kw = ThermalConductivityWater(obj, temperature, sub)
            % Validate the inputs
            arguments
                obj
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            Kw = obj.KwA(1,1,sub)*temperature + obj.KwB(1,1,sub);
        end
        
        %% Specific heat (Cp) function of water (in J/kg/K)
        %   @param temperature in K
        function cp = SpecificHeatWater(obj, temperature, sub)
            % Validate the inputs
            arguments
                obj
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            cp = obj.cpA(1,1,sub)*temperature.^2 + obj.cpB(1,1,sub) * temperature + obj.cpC(1,1,sub);
        end

        %% Diffusion function in water
        %   @param temperature in K
        function D = DiffusionWater(obj, varID, temperature, sub)
            % Validate the inputs
            arguments
                obj
                varID {mustBeInteger} = 1
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            D = obj.refDiff(varID,1,sub)*exp(-1991.805.*((1./temperature)-(1/obj.refDiffTemp(varID,1,sub))));
        end
        
        %% Dispersivity function in water
        %   @param temperature in K
        %   @param ux, uy, uz = velocities in x, y, and z in m/s
        function D = DispersionWater(obj, varID, temperature, ux, uy, uz, sub)
            % Validate the inputs
            arguments
                obj
                varID {mustBeInteger} = 1
                temperature {mustBeNumeric} = 298
                ux {mustBeNumeric} = 0
                uy {mustBeNumeric} = 0
                uz {mustBeNumeric} = 0
                sub {mustBePositive} = 1
            end
            umag = sqrt(ux.^2 + uy.^2 + uz.^2);
            D = obj.DiffusionWater(varID, temperature, sub) + umag*obj.alpha(1,1,sub);
        end
        
        %% Effective Dispersivity function in water
        %   @param temperature in K
        %   @param ux, uy, uz = velocities in x, y, and z in m/s
        function D = EffectiveDispersionWater(obj, varID, temperature, ux, uy, uz, sub)
            % Validate the inputs
            arguments
                obj
                varID {mustBeInteger} = 1
                temperature {mustBeNumeric} = 298
                ux {mustBeNumeric} = 0
                uy {mustBeNumeric} = 0
                uz {mustBeNumeric} = 0
                sub {mustBePositive} = 1
            end
            D = obj.DispersionWater(varID, temperature, ux, uy, uz, sub)*obj.bulk_porosity(1,1,sub)^obj.diff_factor(1,1,sub);
        end
        
        %% Averaged Darcy flow coefficient
        %   units: m^3*s/kg   --> with pressure in Pa = [kg/m/s^2] this
        %   would result in velocities in m/s.
        function K = KozenyCarmannDarcyCoeffient(obj, temperature, sub)
            % Validate the inputs
            arguments
                obj
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            vis = obj.ViscosityWater(temperature, sub);
            K = (obj.particle_diameter(1,1,sub)^2*obj.bulk_porosity(1,1,sub)^3/obj.KozenyCarmannConst(1,1,sub)/(1-obj.bulk_porosity(1,1,sub))^2)./vis;
        end


        %% Calculation of temperature time coefficient
        %   (eb*rho*cpw + (1-eb)*rhos*cps)
        function coeff = TempTimeCoeff(obj, pressure, temperature, sub)
            % Validate the inputs
            arguments
                obj
                pressure    {mustBeNumeric} = 101350
                temperature {mustBeNumeric} = 298
                sub {mustBePositive} = 1
            end
            rho = obj.DensityWater(pressure,temperature,sub);
            cpw = obj.SpecificHeatWater(temperature,sub);
            coeff = obj.bulk_porosity(1,1,sub)*rho.*cpw + (1-obj.bulk_porosity(1,1,sub))*obj.bulk_solids_dens(1,1,sub)*obj.bulk_solids_Cp(1,1,sub);
        end

    end

    %% Private methods
    methods (Access = public)

        %% Function for d_coeff in PDE toolbox
        function dmatrix = d_coeff_fun(obj, sub, location, state)
            Nv = 2+size(obj.mobile_spec_idx,1);   % Variables
            Nl = numel(location.x);               % Locations
            dmatrix = zeros(Nv,Nl);
            dmatrix(1,:) = 0;  % pressure
            dmatrix(2,:) = 1;  % temperature
            for i=3:Nv
                dmatrix(i,:) = 1; % concentration 
            end
        end

        %% Function for d_coeff in PDE toolbox
        function cmatrix = c_coeff_fun(obj, sub, location, state)
            Nv = (2+size(obj.mobile_spec_idx,1))*2;   % Variables (2N form)
            Nl = numel(location.x);               % Locations
            cmatrix = zeros(Nv,Nl);
            cmatrix(1,:) = obj.KozenyCarmannDarcyCoeffient(state.u(2,:),sub);  % pressure
            cmatrix(2,:) = cmatrix(1,:);  % pressure

            cmatrix(3,:) = (obj.bulk_porosity(1,1,sub)*obj.ThermalConductivityWater(state.u(2,:),sub) + ...
                (1-obj.bulk_porosity(1,1,sub)*obj.bulk_solids_cond(1,1,sub))) ... 
                ./obj.TempTimeCoeff(state.u(1,:),state.u(2,:),sub);  % temperature
            cmatrix(4,:) = cmatrix(3,:);  % temperature

            j=1;
            for i=5:2:Nv
                cmatrix(i,:) = obj.EffectiveDispersionWater(j,state.u(2,:),...
                        -obj.KozenyCarmannDarcyCoeffient(state.u(2,:),sub).*state.ux(1,:),...
                        -obj.KozenyCarmannDarcyCoeffient(state.u(2,:),sub).*state.uy(1,:),0,sub) * ...
                        obj.mobile_spec_idx(j,1,sub); % concentration 
                cmatrix(i+1,:) = cmatrix(i,:); % concentration 
                j=j+1;
            end
        end

    end
end

