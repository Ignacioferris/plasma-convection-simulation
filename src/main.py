# visualiser
import FVis3 
import numpy as np

class Convection2D:
    def __init__(self):
        """
        define variables
        """
        # Grid parameters
        self.Nx = 300     # Number of horizontal grid points
        self.Ny = 100     # Number of vertical grid points
        self.Lx = 12.0e6  # Horizontal extent in meters (12 Mm)
        self.Ly = 4.0e6   # Vertical extent in meters (4 Mm)
        self.dx = self.Lx / self.Nx
        self.dy = self.Ly / self.Ny

        # Create meshgrid for coordinates
        self.x = np.linspace(0.5 * self.dx, self.Lx - 0.5 * self.dx, self.Nx)
        self.y = np.linspace(0.5 * self.dy, self.Ly - 0.5 * self.dy, self.Ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij') # Indexing='ij' so X has shape (Nx, Ny)

        # Physical constants
        self.gamma = 5.0 / 3.0       # Adiabatic index for ideal monatomic gas
        self.mu = 0.61               # Mean molecular weight
        self.m_u = 1.660539e-27      # Atomic mass unit in kg
        self.kB = 1.380649e-23       # Boltzmann constant in J/K
        self.G = 6.67430e-11         # Gravitational constant in N m^2/kg^2
        self.M_sun = 1.989e30        # Solar mass in kg
        self.R_sun = 6.957e8         # Solar radius in m
        self.g_abs = self.G * self.M_sun / self.R_sun**2 # Gravitational acceleration
        self.gy = -self.g_abs        # Gravity acts in the negative y-direction (downwards)
        self.T_photosphere = 5778.0  # Temperature of photosphere [K]
        self.P_photosphere = 1.8e4   # Pressure of photosphere [Pa]
        self.nabla = 8             # This can be changed at will


        # Fluid variables (primary variables)
        # These will be Nx x Ny arrays
        self.rho = np.zeros((self.Nx, self.Ny)) # Density (rho) [kg/m^3]
        self.u = np.zeros((self.Nx, self.Ny))   # Horizontal velocity (u) [m/s]
        self.w = np.zeros((self.Nx, self.Ny))   # Vertical velocity (w) [m/s]
        self.e = np.zeros((self.Nx, self.Ny))   # Internal energy density (e) [J/m^3]

        # Initialize arrays for the time derivatives
        self.drho_dt = np.zeros((self.Nx, self.Ny))
        self.du_dt = np.zeros((self.Nx, self.Ny))
        self.dw_dt = np.zeros((self.Nx, self.Ny))
        self.de_dt = np.zeros((self.Nx, self.Ny))
        self.drhou_dt = np.zeros((self.Nx, self.Ny))
        
        # Derived variables
        self.P = np.zeros((self.Nx, self.Ny))   # Pressure (P) [Pa]
        self.T = np.zeros((self.Nx, self.Ny))   # Temperature (T) [K]

        # Time parameters
        self.t = 0.0      # Current time
        self.dt = 0.1     # Timestep
        self.cfl_p = 0.05  # CFL condition factor p

        # Precompute frequently used values
        self.inv_dx = 1.0 / self.dx
        self.inv_dy = 1.0 / self.dy
        self.inv_2dx = 0.5 * self.inv_dx
        self.inv_2dy = 0.5 * self.inv_dy
        self.gamma_minus_1 = self.gamma - 1.0
        self.mu_m_u = self.mu * self.m_u
        self.mu_m_u_g = self.mu_m_u * self.g_abs
        self.kB_inv = 1.0 / self.kB
        self.epsilon = 1e-9

        # A variable to determine if the current data is deleted or not
        self.delete_curr_data = False

        # Gaussian perturbation parameters
        self.apply_perturbation = False # Set to True to add temperature perturbation
        self.apply_isobaric_perturbation = False # Set to True to add temperature perturbation that tries to be isobaric
        self.pert_amplitude_multiplier = 0.8

    def initialise(self):
        """
        initialise temperature, pressure, density and internal energy
        """
        # Initialise P and T based on hydrostatic equilibrium and self.nabla_val
        self.T[:, self.Ny-1] = self.T_photosphere
        self.P[:, self.Ny-1] = self.P_photosphere

        depth = self.Lx - self.y
        temp_gradient_factor = self.nabla * (self.mu * self.m_u * self.g_abs) / self.kB

        # 1 dimensional profile of T and P in terms of the depth
        T_profile_1D = self.T_photosphere + temp_gradient_factor * depth
        P_profile_1D = self.P_photosphere * (T_profile_1D / self.T_photosphere)**(1.0/self.nabla)

        # Fill 2D arrays (assuming homogeneity in x for initial state)
        self.T = np.tile(T_profile_1D[np.newaxis, :], (self.Nx, 1))
        self.P = np.tile(P_profile_1D[np.newaxis, :], (self.Nx, 1))

        # Calculate density and internal energy from P and T
        self.rho = (self.P * self.mu * self.m_u) / (self.kB * self.T)
        self.e = self.P / (self.gamma - 1.0) # Using P = (gamma-1)e 

        # Initial velocity is zero everywhere 
        self.u = np.zeros((self.Nx, self.Ny))
        self.w = np.zeros((self.Nx, self.Ny))

        # Add 2D Gaussian perturbation to temperature if flag is set
        if self.apply_perturbation:
            pert_amplitude = self.pert_amplitude_multiplier * self.T[self.Nx//2, self.Ny//2] # 30% of central temperature
            pert_sigma_x = self.Lx / 5.0
            pert_sigma_y = self.Ly / 5.0
            pert_x0 = self.Lx / 2.0
            pert_y0 = self.Ly / 3.0
            
            gaussian = pert_amplitude * np.exp(
                -((self.X - pert_x0)**2 / (2 * pert_sigma_x**2) +
                  (self.Y - pert_y0)**2 / (2 * pert_sigma_y**2))
            )
            self.T += gaussian
            # Recompute e and P based on perturbed T and original rho (or re-equilibrate P and then e). For consistency, if T changes, e and P should also change.
            self.P = (self.rho * self.kB * self.T) / (self.mu * self.m_u)
            self.e = self.P / (self.gamma - 1.0)

            print("Temperature perturbation applied.")

        # Add 2D Gaussian perturbation to temperature if flag is set
        if self.apply_isobaric_perturbation:
            # Store the original pressure field from hydrostatic equilibrium
            P_hydrostatic = self.P.copy()

            pert_amplitude = self.pert_amplitude_multiplier * self.T[self.Nx//2, self.Ny//2]
            pert_sigma_x = self.Lx / 10.0
            pert_sigma_y = self.Ly / 10.0
            pert_x0 = self.Lx / 2.0
            pert_y0 = self.Ly / 4.0 # Centered lower in the box to give it space to rise
            
            gaussian = pert_amplitude * np.exp(
                -((self.X - pert_x0)**2 / (2 * pert_sigma_x**2) +
                  (self.Y - pert_y0)**2 / (2 * pert_sigma_y**2))
            )
            self.T += gaussian  # Apply temperature perturbation

            # Adjust density to maintain pressure close to the original P_hydrostatic.
            # This makes the heated region less dense and thus buoyant.
            self.rho = (P_hydrostatic * self.mu * self.m_u) / (self.kB * self.T)

            # Recalculate P based on the new T and new rho.
            # This P should now be very close to P_hydrostatic, minimizing pressure waves.
            self.P = (self.rho * self.kB * self.T) / (self.mu * self.m_u)
            
            # Recalculate internal energy density
            self.e = self.P / (self.gamma - 1.0)

            print("Temperature perturbation applied (with isobaric-like density adjustment).")


    def central_x(self, var):
        """
        central difference scheme in x-direction, assuming periodic boundaries
        var_jp1 - var_jm1 / (2*dx)
        """
        # Equivalent to (var[i+1,j] - var[i-1,j]) / (2*dx)
        return (np.roll(var, -1, axis=0) - np.roll(var, 1, axis=0)) * self.inv_2dx

    def central_y(self, var):
        """
        central difference scheme in y-direction
        """
        res = np.zeros_like(var)
    
        # Interior points - standard central difference
        res[:, 1:-1] = (var[:, 2:] - var[:, :-2]) * self.inv_2dy

        # Boundaries - use one-sided 2nd order accurate differences
        # Bottom: forward difference
        res[:, 0] = (-3*var[:, 0] + 4*var[:, 1] - var[:, 2]) * self.inv_2dy
        # Top: backward difference  
        res[:, -1] = (3*var[:, -1] - 4*var[:, -2] + var[:, -3]) * self.inv_2dy
    
        return res

    def upwind_x(self, var, vel_u):
        """
        upwind difference scheme in x-direction
        var - var_im1 / dx if u_ij >=0
        var_ip1 - var / dx if u_ij < 0
        """
        deriv = np.zeros_like(var)
        
        # Periodic boundary conditions in the x-direction are implicitely handled by np.roll
        var_im1 = np.roll(var, 1, axis=0)
        var_ip1 = np.roll(var, -1, axis=0)
        
        deriv[vel_u >= 0] = (var[vel_u >= 0] - var_im1[vel_u >= 0]) * self.inv_dx
        deriv[vel_u < 0]  = (var_ip1[vel_u < 0] - var[vel_u < 0]) * self.inv_dx
        return deriv

    def upwind_y(self, var, vel_w):
        """
        upwind difference scheme in y-direction with proper boundary handling
        """
        deriv = np.zeros_like(var)
    
        # For interior points (j=1 to Ny-2)
        for j in range(1, self.Ny-1):
            if np.any(vel_w[:, j] >= 0):
                # Use backward difference: (var[j] - var[j-1]) / dy
                mask_pos = vel_w[:, j] >= 0
                deriv[mask_pos, j] = (var[mask_pos, j] - var[mask_pos, j-1]) * self.inv_dy
        
            if np.any(vel_w[:, j] < 0):
                # Use forward difference: (var[j+1] - var[j]) / dy
                mask_neg = vel_w[:, j] < 0
                deriv[mask_neg, j] = (var[mask_neg, j+1] - var[mask_neg, j]) * self.inv_dy
    
        # Boundary conditions - use one-sided differences
        # Bottom boundary (j=0)
        deriv[:, 0] = (var[:, 1] - var[:, 0]) * self.inv_dy
    
        # Top boundary (j=Ny-1)
        deriv[:, -1] = (var[:, -1] - var[:, -2]) * self.inv_dy
    
        return deriv

    def boundary_conditions(self):
        """
        boundary conditions for energy, density and velocity
        """
        # Horizontal boundaries (x-direction) are periodic and handled by np.roll in central_x, upwind_x if used directly on self.rho etc.
        
        # Vertical velocity (w): Zero at upper (j=Ny-1) and lower (j=0) boundaries
        self.w[:, 0] = 0.0
        self.w[:, -1] = 0.0

        # The vertical gradient of the horizontal component of the velocity should be zero at the boundary. Solve equations (32) and (33) for du/dy = 0
        self.u[:,  0] = (4*self.u[:, 1] - self.u[:, 2]) / 3.0
        self.u[:, -1] = (4*self.u[:, -2] - self.u[:, -3]) / 3.0

        # Bottom boundary (j=0):
        self.T[:, 0] = self.T[:, 1] # We assume zero gradient for T at the boundaries. I don't know if this is the best wat to approach this.
        self.P[:, 0] = (self.P[:, 1] + 0.5 * self.dy * self.g_abs * self.rho[:, 1]) / (1.0 - 0.5*(self.mu * self.m_u * self.g_abs * self.dy) / (self.kB * self.T[:, 0]))
        self.rho[:, 0] = (self.P[:, 0] * self.mu * self.m_u) / (self.kB * self.T[:, 0])
        self.e[:, 0]   = self.P[:, 0] / (self.gamma - 1.0)

        # Top boundary (j=Ny-1):
        self.T[:, self.Ny-1] = self.T[:, self.Ny-2] # We assume zero gradient for T at the boundaries. I don't know if this is the best wat to approach this.
        self.P[:, self.Ny-1] =  (self.P[:, self.Ny-2] - 0.5 * self.dy * self.g_abs * self.rho[:, self.Ny-2]) / (1.0 + 0.5*(self.mu * self.m_u * self.g_abs * self.dy) / (self.kB * self.T[:, self.Ny-1]))
        self.rho[:, self.Ny-1] = (self.P[:, self.Ny-1] * self.mu * self.m_u) / (self.kB * self.T[:, self.Ny-1])
        self.e[:, self.Ny-1]   = self.P[:, self.Ny-1] / (self.gamma - 1.0)

    # OLD TIMESTEP METHOD
    # def timestep(self):
    #     """
    #     calculate timestep based on CFL condition
    #     """
    #     # Based on equations (25), (26), (27)
    #     epsilon = 1e-9 # Small number to avoid division by zero
    #     velocity_threshold = 1e-5 # Velocities smaller than this are excluded from the CFL condition

    #     # We need to calculate the missing derivatives du/dt and dw/dt, we do so by applying the chain rule
    #     self.du_dt = (self.drhou_dt - self.u * self.drho_dt) / self.rho
    #     self.dw_dt = (self.drhow_dt - self.w * self.drho_dt) / self.rho
        
    #     rel_rho = np.abs(self.drho_dt / (self.rho + epsilon))
    #     rel_u = np.abs(self.du_dt / self.u, where=np.abs(self.u) > velocity_threshold, out=np.zeros_like(self.du_dt))
    #     rel_w = np.abs(self.dw_dt / self.w, where=np.abs(self.w) > velocity_threshold, out=np.zeros_like(self.dw_dt))
    #     rel_e = np.abs(self.de_dt / (self.e + epsilon))
    #     rel_x = np.abs(self.u / self.dx)
    #     rel_y = np.abs(self.w / self.dy)
        
    #     # Find maximum over the grid for each
    #     max_rel_rho = np.max(rel_rho)
    #     max_rel_u = np.max(rel_u)
    #     max_rel_w = np.max(rel_w)
    #     max_rel_e = np.max(rel_e)
    #     max_rel_x = np.max(rel_x)
    #     max_rel_y = np.max(rel_y)

    #     # Overall maximum delta
    #     delta_max = np.max([max_rel_rho, max_rel_u, max_rel_w, max_rel_e, max_rel_x, max_rel_y])
                
    #     # Compute dt
    #     if delta_max == 0:
    #         self.dt = self.cfl_p / (1e-5 + delta_max)  # Default for static case
    #     else:
    #         self.dt = self.cfl_p / delta_max
        
    #     # return 1
    #     print("dt = ", self.dt)
    #     return self.dt

    def timestep(self):
        """
        Calculate timestep based on CFL condition using sound speed and velocities.
        """
        # Compute sound speed
        c = np.sqrt(self.gamma * self.P / self.rho)

        # Compute maximum effective speeds in x and y directions
        max_speed_x = np.max(np.abs(self.u) + c)
        max_speed_y = np.max(np.abs(self.w) + c)
    
        # Compute dt based on CFL condition
        dt_x = self.dx / max_speed_x
        dt_y = self.dy / max_speed_y
        self.dt = self.cfl_p * min(dt_x, dt_y)
    
        # print("dt = ", self.dt)
        return self.dt

    def hydro_solver(self):
        """
        hydrodynamic equations solver
        This function will calculate the time derivatives, then the timestep,
        and finally update the variables.
        """
        du_dx_cen = self.central_x(self.u)
        dw_dy_cen = self.central_y(self.w) # Used in continuity and potentially energy
        drho_dx_up = self.upwind_x(self.rho, self.u)
        drho_dy_up = self.upwind_y(self.rho, self.w)
        self.drho_dt = -self.rho * (du_dx_cen + dw_dy_cen) - self.u * drho_dx_up - self.w * drho_dy_up

        du_dx_up = self.upwind_x(self.u, self.u) # Used in horizontal and vertical momentum
        dw_dy_up = self.upwind_y(self.w, self.w) # Upwind derivative of w w.r.t y
        drhou_dx_up = self.upwind_x(self.rho * self.u, self.u)
        drhou_dy_up = self.upwind_y(self.rho * self.u, self.w)
        dP_dx_cen = self.central_x(self.P)
        self.drhou_dt = - (self.rho * self.u) * (du_dx_up + dw_dy_cen) - self.u * drhou_dx_up - self.w * drhou_dy_up - dP_dx_cen

        # For vertical momentum (see recommendation #2 below for this line)
        dw_dy_up = self.upwind_y(self.w, self.w) # This is (dw/dy)_upwind
        drhow_dy_up = self.upwind_y(self.rho * self.w, self.w)
        drhow_dx_up = self.upwind_x(self.rho * self.w, self.u)
        dP_dy_cen = self.central_y(self.P)
        # Ensure this line for drhow_dt is corrected as per point 2
        # Corrected drhow_dt:
        self.drhow_dt = - (self.rho * self.w) * (dw_dy_up + du_dx_cen) \
                        - self.w * drhow_dy_up \
                        - self.u * drhow_dx_up \
                        - dP_dy_cen + self.rho * self.gy

        deu_dx_up = self.upwind_x(self.e * self.u, self.u)
        dew_dy_up = self.upwind_y(self.e * self.w, self.w)
        de_dx_up = self.upwind_x(self.e, self.u) # Upwind derivative of e w.r.t x
        de_dy_up = self.upwind_y(self.e, self.w) # Upwind derivative of e w.r.t y
        
        self.de_dt = - (self.u * de_dx_up + self.w * de_dy_up) \
                       - (self.e + self.P) * (du_dx_cen + dw_dy_cen)
        
        # --- Timestep Calculation ---
        self.dt = self.timestep() 
        current_dt = self.dt # Store for advancing time self.t
        
        # --- Variable Updates ---
        # Store values at time n that are needed for the RHS of updates if they are modified before use
        rho_n_val = self.rho.copy()  # rho at time n
        u_n_val = self.u.copy()    # u at time n
        w_n_val = self.w.copy()    # w at time n
        e_n_val = self.e.copy()    # e at time n

        # Calculate momentum terms at time n using values at time n
        rhou_n_val = rho_n_val * u_n_val
        rhow_n_val = rho_n_val * w_n_val
        
        # Update density to rho^{n+1}
        self.rho[:] = rho_n_val + self.drho_dt * self.dt
        
        # Update velocities to u^{n+1}, w^{n+1}
        # Denominator self.rho is now rho^{n+1}, which is correct
        self.u[:] = (rhou_n_val + self.drhou_dt * self.dt) / self.rho
        self.w[:] = (rhow_n_val + self.drhow_dt * self.dt) / self.rho
        
        # Update internal energy to e^{n+1}
        self.e[:] = e_n_val + self.de_dt * self.dt

        # Apply boundary conditions to newly updated rho, u, w, e
        self.boundary_conditions()

        # Update derived variables (P, T) from new e and rho
        self.P[:] = (self.gamma - 1.0) * self.e
        self.T[:] = (self.P * self.mu * self.m_u) / (self.kB * self.rho) # Ensure self.rho here is rho^{n+1} and not zero

        # Update current time
        self.t += current_dt

        return current_dt

# --- Main function ---

if __name__ == '__main__':
    sim = Convection2D()
    sim.initialise()

    vis = FVis3.FluidVisualiser(printInfo=True) # Initialize the visualiser

    # --- Parameters for FVis3 ---
    total_sim_time_to_save = 60.0        # Seconds of simulation time to be saved
    simulation_fps_to_save = 0.5          # Save data these times per simulation second

    # Define simulation parameters to save in info.txt
    sim_params_for_fvis = {
        'Nx': sim.Nx, 'Ny': sim.Ny, 
        'dx': sim.dx, 'dy': sim.dy,
        'Lx': sim.Lx, 'Ly': sim.Ly,
        'gamma': sim.gamma, 'mu': sim.mu,
        'g_abs': sim.g_abs,
        'cfl_p': sim.cfl_p,
        'self.nabla_init': 2.0/5.0 + 0.1 # Example if you used this in initialise
    }

    # --- 1. Save Simulation Data ---
    # The arrays (sim.rho, sim.u etc.) are passed by reference and will be updated by hydro_solver
    print(f"Starting data saving for {total_sim_time_to_save} simulation seconds...")
    vis.save_data(
        total_sim_time_to_save, # total simulation time
        sim.hydro_solver,       # update function
        rho=sim.rho,            # Pass by reference
        u=sim.u, 
        w=sim.w, 
        e=sim.e, 
        P=sim.P, 
        T=sim.T,
        sim_fps=simulation_fps_to_save,
        useDblPrec=False,       # Adjust to use double precision, slower
        sim_params=sim_params_for_fvis,
        folder="auto"           # Specify the folder
    )
    print(f"Data saving completed.")

    # --- 2. Animate the Saved Data ---
    # Animate Temperature with velocity quivers
    print(f"Starting animation for Temperature...")
    vis.animate_2D('T',                 # Variable to plot (e.g., 'T', 'rho', 'w', 'v')
        matrixLike=False,               # Set to False if your arrays are (Nx, Ny) and imshow needs (Ny, Nx) after transpose
        folder="default",               # Folder where data was saved
        extent=[0, sim.Lx, 0, sim.Ly],  # Physical domain: [x_min, x_max, y_min, y_max]
        anim_fps=simulation_fps_to_save,# Same fps as the simulation saving
        showQuiver=True,                # Show velocity vectors
        quiverscale=1,                  # Adjust for quiver arrow size
        video_fps=30,                   # Fps of the mp4 video
        N_arrows=20,                    # Number of arrows in quiver plot
        cmap='viridis',                 # Colormap (e.g., 'jet', 'viridis', 'inferno')
        title='2D Convection: Temperature',
        save=True,                      # Save the animation as a video
        video_name='vid_T', 
        units={'T': 'K', 'Lx': 'Mm', 'Lz': 'Mm', 't': 's'} # Lz for vertical axis in FVis3
    )
    print("Temperature animation saved.")

    # Animate Vertical Velocity
    print(f"Starting animation for Vertical Velocity...")
    vis.animate_2D(
        quantity='w',
        matrixLike=False,
        folder="default",
        extent=[0, sim.Lx, 0, sim.Ly],
        showQuiver=False, # Maybe you don't want quivers on a velocity plot? 
        cmap='RdBu_r',    # A diverging colormap good for velocities
        title='2D Convection: Vertical Velocity (w)',
        save=True,
        video_name='vid_w',
        units={'w': 'm/s', 'Lx': 'Mm', 'Lz': 'Mm', 't': 's'}
    )
    print("Vertical velocity animation saved.")

    print(f"Starting animation for Pressure from folder...")
    vis.animate_2D('P',                 # Variable to plot (e.g., 'T', 'rho', 'w', 'v')
        matrixLike=False,               # Set to False if your arrays are (Nx, Ny) and imshow needs (Ny, Nx) after transpose
        folder="default",               # Folder where data was saved
        extent=[0, sim.Lx, 0, sim.Ly],  # Physical domain: [x_min, x_max, y_min, y_max]
        anim_fps=simulation_fps_to_save,# Same fps as the simulation saving
        showQuiver=True,                # Show velocity vectors
        quiverscale=1,                  # Adjust for quiver arrow size
        video_fps=30,                   # Fps of the mp4 video
        N_arrows=20,                    # Number of arrows in quiver plot
        cmap='inferno',                 # Colormap (e.g., 'jet', 'viridis', 'inferno')
        title='2D Convection: Temperature',
        save=True,                      # Save the animation as a video
        video_name='vid_pressure', 
        units={'P': 'Pa', 'Lx': 'Mm', 'Lz': 'Mm', 't': 's'} # Lz for vertical axis in FVis3
    )
    print("Pressure animation saved.")

    print(f"Starting animation for Energy density from folder...")
    vis.animate_2D('e',                 # Variable to plot (e.g., 'T', 'rho', 'w', 'v')
        matrixLike=False,               # Set to False if your arrays are (Nx, Ny) and imshow needs (Ny, Nx) after transpose
        folder="default",               # Folder where data was saved
        extent=[0, sim.Lx, 0, sim.Ly],  # Physical domain: [x_min, x_max, y_min, y_max]
        anim_fps=simulation_fps_to_save,# Same fps as the simulation saving
        showQuiver=True,                # Show velocity vectors
        quiverscale=1,                  # Adjust for quiver arrow size
        video_fps=30,                   # Fps of the mp4 video
        N_arrows=20,                    # Number of arrows in quiver plot
        cmap='inferno',                 # Colormap (e.g., 'jet', 'viridis', 'inferno')
        title='2D Convection: Temperature',
        save=True,                      # Save the animation as a video
        video_name='vid_energy_density', 
        units={'e': 'J/m^3', 'Lx': 'Mm', 'Lz': 'Mm', 't': 's'} # Lz for vertical axis in FVis3
    )
    print("Energy density animation saved.")

    if sim.delete_curr_data:
        vis.delete_current_data()