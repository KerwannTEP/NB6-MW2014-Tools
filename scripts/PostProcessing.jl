using DelimitedFiles
using ArgParse
using Glob
using Plots
using LaTeXStrings
using Plots.PlotMeasures
using HDF5
using SpecialFunctions
using Statistics

tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--isnap"
    help = "Snapshot index for IOM analysis."
    arg_type = Int64
    default = 1
    "--name"
    help = "Name of the model. Default: eccentric"
    arg_type = String
    default = "eccentric"
    "--rmax"
    help = "Position bounds (in kpc)."
    arg_type = Float64
    default = 25.0
end

parsed_args = parse_args(tabargs)
const isnap = parsed_args["isnap"]
const name_model = parsed_args["name"]
const rmax = parsed_args["rmax"]

const G_in_pc1_MSun1_kms2 = 4.3009172706e-3

const path_to_data = "path/to/data/" # Location of the folder containing the data
const path_to_plot = "path/to/plot/" # Location of the folder containing the plots

################################################################################################################
# MWPotential2014 (Bovy 2015)
# https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract
################################################################################################################

# Copy MWpotential.f file from NBODY6++GPU

# Bulge
const b_alpha = 1.8
const b_rc = 1.9 # kpc 
const b_norm = 0.05 # Contribution from the bulge to the force at R0=8 kpc

# Disk 
const d_a = 3.0 # kpc 
const d_b = 0.28 # kpc 
const d_norm = 0.6 # Contribution from the disk to the force at R0=8 kpc

# Halo 
const h_a = 16.0 # kpc 
const h_norm = 0.35  # Contribution from the halo to the force at R0=8 kpc

# Scale : Force
const F_scale = 6.32793804994 # pc/Myr^2 : Force from the galaxy at R0=8 kpc


# Forces

# fpowercut.f
# This returns gamma(a) * gammad(x, a) from nbody6++
function gamma_low(x::Float64, a::Float64)

    return gamma(a) - gamma(a, x)

end

function Menc_bulge(r::Float64)

    a = 1.5 - b_alpha/2.0
    xup = (r/b_rc)^2
    ga = gamma_low(xup, a) # = gamma(a) * gammad(x, a)
    b_amp = 1.0
    menc = b_amp*2*pi*b_rc^(3-b_alpha)*ga

    return menc
end

function force_bulge_at_R0()
    x = 8.0 # kpc 
    y = 0.0 # kpc 
    z = 0.0 # kpc 
    r = x
    Fr = - Menc_bulge(r)/r^2
    
    return Fr
end


# miyamoto.f
function force_disk_at_R0()
    x = 8.0 # kpc 
    y = 0.0 # kpc 
    z = 0.0 # kpc 
    R = x
    bz = sqrt(d_b^2 + z^2)
    az = sqrt(R^2 + (d_a+bz)^2)
    d_amp = 1.0
    az3 = d_amp/az^3
    FR = -az3*x

    return FR

end

# fnfw.f
function force_halo_at_R0()
    x = 8.0 # kpc 
    y = 0.0 # kpc 
    z = 0.0 # kpc 
    r = x
    ra = h_a + r
    h_amp = 1.0
    dp = h_amp*(1/(r^2*ra) - log(1+r/h_a)/r^3)
    Fr = dp*x

    return Fr 
end

function normalize_amplitude()

    bamp = b_norm * abs(F_scale/force_bulge_at_R0())
    damp = d_norm * abs(F_scale/force_disk_at_R0())
    hamp = h_norm * abs(F_scale/force_halo_at_R0())

    return bamp, damp, hamp
end

const b_amp, d_amp, h_amp = normalize_amplitude()


# Potentials
# Potential unit: kpc*pc/Myr^2

function potential_bulge(r::Float64, amp::Float64=b_amp)

    xup = (r/b_rc)^2.0
    a = 1.0 - 0.5*b_alpha
    a2 = 1.5 - 0.5*b_alpha
    g1 = gamma_low(xup, a) # gammad(xup, a) * gamma(a)
    g2 = gamma_low(xup, a2) # gammad(xup, a2) * gamma(a2)

    pot = amp*2*pi*b_rc^(3.0-b_alpha)*(g1/b_rc-g2/r)

    return pot 
end

function potential_disk(R::Float64, z::Float64, amp::Float64=d_amp)

    pot = -amp/sqrt(R^2+(d_a + sqrt(z^2+d_b^2))^2)

    return pot

end

function potential_halo(r::Float64, amp::Float64=h_amp)

    pot = -amp*log(1.0 + r/h_a)/r 

    return pot

end

function E_max(amp::Float64=b_amp)

    a = 1.0 - 0.5*b_alpha
    g1 = gamma(a) # gammad(xup, a) * gamma(a)

    pot = amp*2*pi*b_rc^(3.0-b_alpha)*(g1/b_rc)

    return pot 
end


################################################################################################################
# Save important quantities and plots
################################################################################################################

function get_list_files()

    list_files = glob(path_to_data * name_model * "/output/out_*.txt")
    nbt = length(list_files)
    list_time = zeros(Int64, nbt)

    for i=1:nbt 
        time = parse(Int64, split(split(list_files[i], "_")[end],".")[1])
        list_time[i] = time
    end
    p = sortperm(list_time)

    return list_files[p], nbt
end

function treat_data()

    ###############################################
    # Read data
    ###############################################

    data_global = readdlm(path_to_data * name_model * "/global.30")
    data_lag = readdlm(path_to_data * name_model * "/lagr.7")

    list_files, nbt = get_list_files()

    tab_Rc_Vc = zeros(Float64, nbt, 6)
    tab_dist_cluster_in_kpc = zeros(Float64, nbt)
    tab_velocity_cluster_in_kms = zeros(Float64, nbt)

    tab_time_in_HU = zeros(Float64, nbt)
    tab_E_in_HU = zeros(Float64, nbt)
    tab_L_in_HU = zeros(Float64, nbt, 3)
    tab_virial_ratio = zeros(Float64, nbt)

    tab_lag_in_HU = zeros(Float64, nbt, 5)
    tab_frac_unbound = zeros(Float64, nbt)

    # Get units
    data0 = readdlm(list_files[1])
    data0_info = Float64.(data0[1, :])
    data0_star = Float64.(data0[2:end, 1:8])

    pc_per_HU = data0_info[2]
    kpc_per_HU = pc_per_HU/1000.0
    Msun_per_HU = data0_info[3]
    Myr_per_HU = data0_info[8]
    kms_per_HU = data0_info[9]
    Npart = length(data0_star[:, 1])

    Nub = 0
    tab_E_wrt_host_in_kms = zeros(Float64, Npart)
    tab_Lz_wrt_host_in_kpckms = zeros(Float64, Npart)
    tab_E_wrt_host_in_kms_unbound = zeros(Float64, Npart)
    tab_Lz_wrt_host_in_kpckms_unbound = zeros(Float64, Npart)
    tab_E_wrt_cluster_in_kms = zeros(Float64, Npart)

    tab_pos_unbound_in_kpc = zeros(Float64, Npart, 3)
    tab_pos_snapshot_in_kpc = zeros(Float64, Npart, 3)
    tab_vel_snapshot_in_kms = zeros(Float64, Npart, 3)
    tab_phi1_phi2_in_deg = zeros(Float64, Npart, 2)

    tab_rc_norm = zeros(Float64, 3)
    tab_tc_norm = zeros(Float64, 3)
    tab_Lc_norm = zeros(Float64, 3)

    tab_pos_corotating_in_kpc = zeros(Float64, Npart, 3)
    
    # temp 
    xc_snap = 0.0
    yc_snap = 0.0
    zc_snap = 0.0

    Threads.@threads for i=1:nbt
        namefile = list_files[i]
        data = readdlm(namefile)

        data_info = Float64.(data[1, :])
        data_star = Float64.(data[2:end, 1:8])

        # Cluster's position (in HU)
        xcl = data_info[13]
        ycl = data_info[14]
        zcl = data_info[15]
        vxcl = data_info[16]
        vycl = data_info[17]
        vzcl = data_info[18] 

        tab_Rc_Vc[i, 1] = xcl * kpc_per_HU
        tab_Rc_Vc[i, 2] = ycl * kpc_per_HU
        tab_Rc_Vc[i, 3] = zcl * kpc_per_HU
        tab_Rc_Vc[i, 4] = vxcl * kms_per_HU
        tab_Rc_Vc[i, 5] = vycl * kms_per_HU
        tab_Rc_Vc[i, 6] = vzcl * kms_per_HU

        tab_dist_cluster_in_kpc[i] = sqrt(xcl^2 + ycl^2 + zcl^2) .* kpc_per_HU
        tab_velocity_cluster_in_kms[i] = sqrt(vxcl^2 + vycl^2 + vzcl^2) .* kms_per_HU

        # Read stellar data
        tab_time_in_HU[i] = data_info[1]
        r_dens_in_HU = data_info[5:7]

        # Read cluster information
        tab_m_in_HU = data_star[:, 1]
        tab_r_in_HU = data_star[:, 2:4]
        tab_v_in_HU = data_star[:, 5:7]
        tab_pot_in_HU = data_star[:, 8]

        Untot = 0.0 # Nbody
        Uhtot = 0.0 # Host
        Ktot = 0.0

        Ktot_c = 0.0

        Lxtot = 0.0
        Lytot = 0.0
        Lztot = 0.0

        for k=1:Npart 
            xc, yc, zc = tab_r_in_HU[k, :]
            x = xc + xcl
            y = yc + ycl
            z = zc + zcl
            R = sqrt(x^2 + y^2)
            r = sqrt(R^2 + z^2)

            vxc, vyc, vzc = tab_v_in_HU[k, :]
            vx = vxc + vxcl
            vy = vyc + vycl
            vz = vzc + vzcl

            m = tab_m_in_HU[k]
            phin = tab_pot_in_HU[k]

            # Potential unit: kpc*pc/Myr^2
            # Multiply by 1/kpc_per_HU 1/pc_per_HU Myr_per_HU^2 to get HU
            factor_pot = (1.0/kpc_per_HU * 1.0/pc_per_HU * Myr_per_HU^2)

            # Bulge
            psib = potential_bulge(r*kpc_per_HU) * factor_pot

            # Disk
            psid = potential_disk(R*kpc_per_HU, z*kpc_per_HU) * factor_pot

            # Halo
            psih = potential_halo(r*kpc_per_HU) * factor_pot

            
            # Total galactic potential 
            psi_gal = psib + psid + psih

            # Potential energy of particle
            Uni = m * phin
            Uhi = m * psi_gal
            Untot += Uni 
            Uhtot += Uhi 

            # Kinetic energy of particle
            Ki = 0.5 * m * (vx^2 + vy^2 + vz^2)
            Ktot += Ki 

            # Angular momentum
            Lxi = m*(y*vz-z*vy)
            Lyi = m*(z*vx-x*vz)
            Lzi = m*(x*vy-y*vx) 

            Lxtot += Lxi
            Lytot += Lyi
            Lztot += Lzi

            # Energy within the cluster frame 
            Ui_c = Uni
            Ki_c = 0.5 * m * (vxc^2 + vyc^2 + vzc^2)
            Ei_c = Ki_c + Ui_c 
            Ktot_c += Ki_c
            
            # Unbound particles at snapshot isnap
            if (i == isnap)
                if (Ei_c >= 0.0)
                    Nub += 1
                end
            end


        end

        Utot = 0.5 * Untot + Uhtot
        Etot = Ktot + Utot
        tab_E_in_HU[i] = Etot 

        tab_virial_ratio[i] = 2*Ktot_c/abs(0.5 * Untot)

        tab_L_in_HU[i, 1] = Lxtot
        tab_L_in_HU[i, 2] = Lytot
        tab_L_in_HU[i, 3] = Lztot


        # Snapshot IOM
        # Use specific IOMs (i.e. IOMs per unit mass)
        # Compute a map of (phi1, phi2)
        if (i == isnap)

            # Compute co-rotating frame (rcl, tcl, Lcl) -> (xs,ys,zs)
            # Lcl = rcl x vcl 
            # tcl = Lcl x rcl

            # Lcl_x = ycl vzcl - zcl vycl
            # Lcl_y = zcl vxcl - xcl vzcl
            # Lcl_z = xcl vycl - ycl vxcl

            # tcl_x = Lcl_y zcl - Lcl_z ycl
            # tcl_y = Lcl_z xcl - Lcl_x zcl
            # tcl_z = Lcl_x ycl - Lcl_y xcl

            # r = x  e_x   + y  e_y   + z  e_z
            # r = x' t{rcl} + y' t{tcl} + z' t{Lcl}

            # xs = <r,rcl>/|rcl| = (x xcl + y ycl + z zcl)/|rcl|
            # ys = <r,tcl>/|tcl| = (x tcl_x + y tcl_y + z tcl_z)/|tcl|
            # zs = <r,Lcl>/|Lcl| = (x Lcl_x + y Lcl_y + z Lcl_z)/|Lcl|

            xc_snap = xcl * kpc_per_HU
            yc_snap = ycl * kpc_per_HU
            zc_snap = zcl * kpc_per_HU

            rcl = sqrt(xcl^2 + ycl^2 + zcl^2)

            Lcl_x = ycl * vzcl - zcl * vycl
            Lcl_y = zcl * vxcl - xcl * vzcl
            Lcl_z = xcl * vycl - ycl * vxcl
            Lcl = sqrt(Lcl_x^2 + Lcl_y^2 + Lcl_z^2)

            tcl_x = Lcl_y * zcl - Lcl_z * ycl
            tcl_y = Lcl_z * xcl - Lcl_x * zcl
            tcl_z = Lcl_x * ycl - Lcl_y * xcl
            tcl = sqrt(tcl_x^2 + tcl_y^2 + tcl_z^2)

            tab_rc_norm[1] = xcl/rcl * 5
            tab_rc_norm[2] = ycl/rcl * 5
            tab_rc_norm[3] = zcl/rcl * 5

            tab_tc_norm[1] = tcl_x/tcl * 5
            tab_tc_norm[2] = tcl_y/tcl * 5
            tab_tc_norm[3] = tcl_z/tcl * 5

            tab_Lc_norm[1] = Lcl_x/Lcl * 5
            tab_Lc_norm[2] = Lcl_y/Lcl * 5
            tab_Lc_norm[3] = Lcl_z/Lcl * 5

            index = 1
            for k=1:Npart 
                xc, yc, zc = tab_r_in_HU[k, :]
                x = xc + xcl
                y = yc + ycl
                z = zc + zcl
                R = sqrt(x^2 + y^2)
                r = sqrt(R^2 + z^2)

                vxc, vyc, vzc = tab_v_in_HU[k, :]
                vx = vxc + vxcl
                vy = vyc + vycl
                vz = vzc + vzcl

                phin = tab_pot_in_HU[k]

                # Potential unit: kpc*pc/Myr^2
                # Multiply by 1/kpc_per_HU 1/pc_per_HU Myr_per_HU^2 to get HU
                factor_pot = (1.0/kpc_per_HU * 1.0/pc_per_HU * Myr_per_HU^2)

                # Bulge
                psib = potential_bulge(r*kpc_per_HU) * factor_pot

                # Disk
                psid = potential_disk(R*kpc_per_HU, z*kpc_per_HU) * factor_pot

                # Halo
                psih = potential_halo(r*kpc_per_HU) * factor_pot


                # Total galactic potential 
                psi_gal = psib + psid + psih

                # Potential energy of particle
                Ui = (phin + psi_gal)

                # Kinetic energy of particle
                Ki = 0.5 * (vx^2 + vy^2 + vz^2)

                # (Specific) Angular momentum
                Lzi_in_HU = (x*vy-y*vx) 
                
                Ui_c = phin
                Ki_c = 0.5 * (vxc^2 + vyc^2 + vzc^2)
                Ei_c = Ki_c + Ui_c 

                Ei_in_HU = Ki + Ui

                # Coordinates in the sky place 
                xs = (x * xcl   + y * ycl   + z * zcl)/rcl
                ys = (x * tcl_x + y * tcl_y + z * tcl_z)/tcl
                zs = (x * Lcl_x + y * Lcl_y + z * Lcl_z)/Lcl

                tab_pos_corotating_in_kpc[k, 1] = xs * kpc_per_HU
                tab_pos_corotating_in_kpc[k, 2] = ys * kpc_per_HU
                tab_pos_corotating_in_kpc[k, 3] = zs * kpc_per_HU

                phi1 = atan(ys, xs)
                phi2 = asin(zs/sqrt(xs^2+ys^2+zs^2))

                tab_phi1_phi2_in_deg[k, 1] = phi1 * 180.0/pi
                tab_phi1_phi2_in_deg[k, 2] = phi2 * 180.0/pi

                tab_pos_snapshot_in_kpc[k, 1] = x * kpc_per_HU
                tab_pos_snapshot_in_kpc[k, 2] = y * kpc_per_HU
                tab_pos_snapshot_in_kpc[k, 3] = z * kpc_per_HU
                tab_vel_snapshot_in_kms[k, 1] = vx * kms_per_HU
                tab_vel_snapshot_in_kms[k, 2] = vy * kms_per_HU
                tab_vel_snapshot_in_kms[k, 3] = vz * kms_per_HU

                tab_E_wrt_cluster_in_kms[k] = Ei_c * (kms_per_HU)^2

                tab_E_wrt_host_in_kms[k] = Ei_in_HU * (kms_per_HU)^2
                tab_Lz_wrt_host_in_kpckms[k] = Lzi_in_HU * kpc_per_HU * kms_per_HU


                if (Ei_c >= 0.0)

                    # Convert to astrophysical units 
                    tab_E_wrt_host_in_kms_unbound[index] = Ei_in_HU * (kms_per_HU)^2
                    tab_Lz_wrt_host_in_kpckms_unbound[index] = Lzi_in_HU * kpc_per_HU * kms_per_HU
                    tab_pos_unbound_in_kpc[index, 1] = x * kpc_per_HU
                    tab_pos_unbound_in_kpc[index, 2] = y * kpc_per_HU
                    tab_pos_unbound_in_kpc[index, 3] = z * kpc_per_HU

                    index += 1
                end
                

            end

        end

    end

    println("==============================")

    # Energy of the cluster's COM
    factor_pot = (1.0/kpc_per_HU * 1.0/pc_per_HU * Myr_per_HU^2)

    Kb = 0.5 * tab_velocity_cluster_in_kms[1]^2 # in (km/s)^2
    Ub = (potential_bulge(tab_dist_cluster_in_kpc[1]) 
        + potential_disk(sqrt(tab_Rc_Vc[1, 1]^2+tab_Rc_Vc[1, 2]^2), tab_Rc_Vc[1, 3])
        + potential_halo(tab_dist_cluster_in_kpc[1]))  * factor_pot * kms_per_HU^2  # in (km/s)^2

    Eb = Kb + Ub

    println("Barycentric energy of the cluster : ", Eb, " (km/2)^2")
    
    E_esc = E_max() * factor_pot * kms_per_HU^2  # in (km/s)^2
    println("Maximum bound energy              : ", E_esc, " (km/2)^2")

    if (Eb < E_esc)
        println("-> The cluster is bound.")
    else
        println("-> The cluster is unbound.")
    end

    println("==============================")

    ###############################################
    # Save in HDF5 file
    ###############################################

    time_snapshot_in_Myr = tab_time_in_HU[isnap] .* Myr_per_HU
    time_snapshot_in_Myr = round(time_snapshot_in_Myr, digits=1) # Cutoff digits

    # Long-term relaxation

    mkpath(path_to_data * "post_treatment/"*name_model*"/")
    namefile_hf5 = path_to_data * "post_treatment/"*name_model*"/iom_cluster_"*name_model*".hf5"
    file = h5open(namefile_hf5, "w")

    write(file, "data_Etot_wrt_host", tab_E_in_HU)
    write(file, "data_Rc_kpc", tab_Rc_Vc[:, 1:3])
    write(file, "data_Vc_kms", tab_Rc_Vc[:, 4:6])
    write(file, "data_virial_ratio_cluster", tab_virial_ratio)

    write(file, "data_time_HU", tab_time_in_HU)
    write(file, "data_time_Myr", tab_time_in_HU .* Myr_per_HU)
    
    write(file, "data_Lx_wrt_host", tab_L_in_HU[:, 1])
    write(file, "data_Ly_wrt_host", tab_L_in_HU[:, 2])
    write(file, "data_Lz_wrt_host", tab_L_in_HU[:, 3])

    write(file, "data_unbound_frac", data_global[2:end,31] ./ Npart)

    write(file, "data_lagrange_rad_00.1_pc", data_lag[3:end, 2] .* pc_per_HU)
    write(file, "data_lagrange_rad_00.3_pc", data_lag[3:end, 3] .* pc_per_HU)
    write(file, "data_lagrange_rad_00.5_pc", data_lag[3:end, 4] .* pc_per_HU)
    write(file, "data_lagrange_rad_01_pc", data_lag[3:end, 5] .* pc_per_HU)
    write(file, "data_lagrange_rad_03_pc", data_lag[3:end, 6] .* pc_per_HU)
    write(file, "data_lagrange_rad_05_pc", data_lag[3:end, 7] .* pc_per_HU)
    write(file, "data_lagrange_rad_10_pc", data_lag[3:end, 8] .* pc_per_HU)
    write(file, "data_lagrange_rad_20_pc", data_lag[3:end, 9] .* pc_per_HU)
    write(file, "data_lagrange_rad_30_pc", data_lag[3:end, 10] .* pc_per_HU)
    write(file, "data_lagrange_rad_40_pc", data_lag[3:end, 11] .* pc_per_HU)
    write(file, "data_lagrange_rad_50_pc", data_lag[3:end, 12] .* pc_per_HU)
    write(file, "data_lagrange_rad_60_pc", data_lag[3:end, 13] .* pc_per_HU)
    write(file, "data_lagrange_rad_70_pc", data_lag[3:end, 14] .* pc_per_HU)
    write(file, "data_lagrange_rad_80_pc", data_lag[3:end, 15] .* pc_per_HU)
    write(file, "data_lagrange_rad_90_pc", data_lag[3:end, 16] .* pc_per_HU)
    write(file, "data_lagrange_rad_95_pc", data_lag[3:end, 17] .* pc_per_HU)
    write(file, "data_lagrange_rad_99_pc", data_lag[3:end, 18] .* pc_per_HU)
    write(file, "data_lagrange_rad_100_pc", data_lag[3:end, 19] .* pc_per_HU)
    write(file, "data_r_core_pc", data_lag[3:end, 20] .* pc_per_HU)

    write(file, "Npart", Npart)
    write(file, "pc_per_HU", pc_per_HU)
    write(file, "kpc_per_HU", kpc_per_HU)
    write(file, "G_in_pc1_MSun1_kms2", G_in_pc1_MSun1_kms2)
    write(file, "Msun_per_HU", Msun_per_HU)

    write(file, "Myr_per_HU", Myr_per_HU)
    write(file, "kms_per_HU", kms_per_HU)

    close(file)


    # Snapshot

    mkpath(path_to_data * "post_treatment/"*name_model*"/")
    namefile_hf5 = path_to_data * "post_treatment/"*name_model*"/snapshot_cluster_"*name_model*"_time_"*string(time_snapshot_in_Myr)*"_Myr.hf5"
    file = h5open(namefile_hf5, "w")

    write(file, "data_phi1_phi2_in_deg",tab_phi1_phi2_in_deg)
    write(file, "data_E_wrt_host_in_kms_unbound", tab_E_wrt_host_in_kms_unbound[1:Nub])
    write(file, "data_Lz_wrt_host_in_kpckms_unbound", tab_Lz_wrt_host_in_kpckms_unbound[1:Nub])

    write(file, "data_E_wrt_cluster_in_kms", tab_E_wrt_cluster_in_kms)
    write(file, "data_E_wrt_host_in_kms", tab_E_wrt_host_in_kms)
    write(file, "data_Lz_wrt_host_in_kpckms", tab_Lz_wrt_host_in_kpckms)

    write(file, "data_pos_in_kpc", tab_pos_snapshot_in_kpc)
    write(file, "data_vel_in_kms", tab_vel_snapshot_in_kms)

    write(file, "Npart", Npart)
    write(file, "pc_per_HU", pc_per_HU)
    write(file, "kpc_per_HU", kpc_per_HU)
    write(file, "G_in_pc1_MSun1_kms2", G_in_pc1_MSun1_kms2)
    write(file, "Msun_per_HU", Msun_per_HU)

    write(file, "Myr_per_HU", Myr_per_HU)
    write(file, "kms_per_HU", kms_per_HU)

    close(file)


    ###############################################
    # Plot data
    ###############################################

    # Total energy
    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_E_in_HU, 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Energy [HU]", 
        # xticks=0:500:5000,
        xminorticks=2,
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/energy.pdf"
    savefig(plt, namefile_pdf)

    # Fractional (total) energy
    dataFracE = abs.(1.0 .- tab_E_in_HU[2:end] ./ tab_E_in_HU[1])
    minE = 10.0^floor(Int64, max(log10(minimum(dataFracE)), -16))
    maxE = 10.0^(floor(Int64, log10(maximum(dataFracE)))+1)

    plt = plot(tab_time_in_HU[2:end] .* Myr_per_HU, dataFracE, 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fractional energy", 
        yaxis=:log10,
        # xticks=0:500:5000,
        xminorticks=2,
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        ylims=(minE, maxE),
        frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/frac_energy.pdf"
    savefig(plt, namefile_pdf)

    # Total angular momentum Lz
    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_L_in_HU[:, 3], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel=L"L_z"*" [HU]", 
        # xticks=0:500:5000,
        xminorticks=2,
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/momentum_z.pdf"
    savefig(plt, namefile_pdf)

    # Fractional angular momentum Lz
    dataFracLz = abs.(1.0 .- tab_L_in_HU[2:end, 3] ./ tab_L_in_HU[1,3])
    minLz = 10.0^floor(Int64, max(log10(minimum(dataFracLz)), -16))
    maxLz = 10.0^(floor(Int64, log10(maximum(dataFracLz)))+1)

    plt = plot(tab_time_in_HU[2:end] .* Myr_per_HU, dataFracLz, 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fractional angular momentum", 
        yaxis=:log10,
        # xticks=0:500:5000,
        xminorticks=2,
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        ylims=(minLz, maxLz),
        frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/frac_momentum_z.pdf"
    savefig(plt, namefile_pdf)

    # Virial ratio
    plt = plot(tab_time_in_HU[1:end] .* Myr_per_HU, tab_virial_ratio, 
        labels=:false, title="Virial ratio",
        xlabel="Time [Myr]", 
        ylabel=L"2 K/|U_{\mathrm{c}}|", 
        yaxis=:log10,
        # xticks=0:500:5000,
        xminorticks=2,
        yticks=10.0 .^ (-20:1:5),
        yminorticks=10,
        xlims=(tab_time_in_HU[1], tab_time_in_HU[end])  .* Myr_per_HU,
        frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/virial_ratio.pdf"
    savefig(plt, namefile_pdf)

    # Distance to cluster
    dmax = maximum(tab_dist_cluster_in_kpc) * 1.1
    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_dist_cluster_in_kpc,
                labels=false,
                xlabel="Time [Myr]", ylabel="Distance to cluster [kpc]",
                xlims=(tab_time_in_HU[1], tab_time_in_HU[end]) .* Myr_per_HU,
                # xticks=0:500:5000,
                xminorticks=2,
                # ylims=(7.9995, 8.0075),
                ylims=(0.0, dmax),
                marker=true, markersize=2,
                frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/distance_cluster.pdf"
    savefig(plt, namefile_pdf)

    # Velocity of the cluster
    dmax = maximum(tab_velocity_cluster_in_kms) * 1.1
    plt = plot(tab_time_in_HU .* Myr_per_HU, tab_velocity_cluster_in_kms,
                labels=false,
                xlabel="Time [Myr]", ylabel="Cluster velocity [km/s]",
                xlims=(tab_time_in_HU[1], tab_time_in_HU[end]) .* Myr_per_HU,
                # xticks=0:500:5000,
                xminorticks=2,
                # ylims=(0.0, dmax),
                marker=true, markersize=2,
                frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/velocity_cluster.pdf"
    savefig(plt, namefile_pdf)

    # Unbound particles
    plt = plot(data_global[2:end,1] .* Myr_per_HU, data_global[2:end,31] ./Npart .* 100.0,
                labels=false,
                xlabel="Time [Myr]", ylabel="Fraction of unbound stars [%]",
                xlims=(data_global[2,1], data_global[end,1]) .* Myr_per_HU,
                # xticks=0:500:5000,
                xminorticks=2,
                # yticks=0:10:100,
                yminorticks=5,
                frame=:box)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/unbound_fraction.pdf"
    savefig(plt, namefile_pdf)

    # Lagrange radii and core radius

    plt = plot(data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,20] .* pc_per_HU, yaxis=:log10, label=L"R_{c}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,5] .* pc_per_HU, yaxis=:log10, label=L"r_{01}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,8] .* pc_per_HU, yaxis=:log10, label=L"r_{10}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,9] .* pc_per_HU, yaxis=:log10, label=L"r_{20}")
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,12] .* pc_per_HU, yaxis=:log10, label=L"r_{50}", legend=:bottomleft)
    plot!(plt, data_lag[3:end,1] .* Myr_per_HU, data_lag[3:end,16] .* pc_per_HU, yaxis=:log10, label=L"r_{90}", legend=:topleft)
    plot!(plt, legend=:topleft, xlabel="Time [Myr]", ylabel="Radii [pc]", frame=:box)
    plot!(plt, yticks=10.0 .^ (-5:1:10), yminorticks=10)
    plot!(plt, xticks=0:500:5000, xminorticks=2)
    plot!(plt, xlims=(data_lag[3,1], data_lag[end,1]) .* Myr_per_HU)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/lagrange_radii.pdf"
    savefig(plt, namefile_pdf)


    # Snapshot IOM
    meanE = mean(tab_E_wrt_host_in_kms_unbound[1:Nub])
    varE = var(tab_E_wrt_host_in_kms_unbound[1:Nub], corrected=true, mean=meanE)
    sigmaE = sqrt(varE)
    tabDeltaE =  (tab_E_wrt_host_in_kms_unbound[1:Nub] .- meanE) ./ sigmaE

    meanLz = mean(tab_Lz_wrt_host_in_kpckms_unbound[1:Nub])
    varLz = var(tab_Lz_wrt_host_in_kpckms_unbound[1:Nub], corrected=true, mean=meanLz)
    sigmaLz = sqrt(varLz)
    tabDeltaLz =  (tab_Lz_wrt_host_in_kpckms_unbound[1:Nub] .- meanLz) ./ sigmaLz

    println("<E>  [unbound, 10^4 (km/s)^2  ] = ", meanE/10^4)
    println("<Lz> [unbound, 10^2 kpc (km/s)] = ", meanLz/10^2)

    

    s = 1.0

    plt = scatter(tab_pos_unbound_in_kpc[1:Nub, 1],
                tab_pos_unbound_in_kpc[1:Nub, 2],
                tab_pos_unbound_in_kpc[1:Nub, 3],
                xlims=(-rmax, rmax), ylims=(-rmax, rmax), zlims=(-rmax, rmax),
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", zlabel=L"z"*" [kpc]", 
                framestyle=:box, labels=false,
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                markerstrokewidth = 0,
                camera=(35, 10),
                title="t = "*string(time_snapshot_in_Myr)*" Myr")

    scatter!(plt, [0], [0], [0], label=false, color=:red)

    plot!(plt, [xc_snap, xc_snap+tab_rc_norm[1]],
               [yc_snap, yc_snap+tab_rc_norm[2]], 
               [zc_snap, zc_snap+tab_rc_norm[3]],
                arrow=true, color=:cyan, linewidth=1,
                label=L"\mathbf{\hat{r}}_{\mathrm{c}}")

    plot!(plt, [xc_snap, xc_snap+tab_tc_norm[1]],
               [yc_snap, yc_snap+tab_tc_norm[2]], 
               [zc_snap, zc_snap+tab_tc_norm[3]],
                arrow=true, color=:yellow, linewidth=1,
                label=L"\mathbf{\hat{t}}_{\mathrm{c}}")

    plot!(plt, [xc_snap, xc_snap+tab_Lc_norm[1]],
               [yc_snap, yc_snap+tab_Lc_norm[2]], 
               [zc_snap, zc_snap+tab_Lc_norm[3]],
               arrow=true, color=:magenta, linewidth=1,
               label=L"\mathbf{\hat{L}}_{\mathrm{c}}")

    display(plt)
    readline()


    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/snapshot_stream_t_"*string(time_snapshot_in_Myr)*"_Myr.pdf"
    savefig(plt, namefile_pdf)

    plt = scatter(tab_pos_corotating_in_kpc[:, 1],
                tab_pos_corotating_in_kpc[:, 2],
                xlims=(-rmax, rmax), ylims=(-rmax, rmax), 
                xlabel=L"x_{\mathrm{corr}}"*" [kpc]", ylabel=L"y_{\mathrm{corr}}"*" [kpc]", 
                framestyle=:box, labels=false,
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                markerstrokewidth = 0,
                title="t = "*string(time_snapshot_in_Myr)*" Myr")

    scatter!(plt, [0], [0], label=false, color=:red)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/snapshot_corotating_stream_t_"*string(time_snapshot_in_Myr)*"_Myr.pdf"
    savefig(plt, namefile_pdf)


    s = 1.0
    # ang_max = 180.0 # deg

    plt = scatter(tab_phi1_phi2_in_deg[:, 1], tab_phi1_phi2_in_deg[:, 2],
                # xlims=(-ang_max, ang_max), ylims=(-ang_max, ang_max)
                xlabel=L"\phi_1"*" [deg]", ylabel=L"\phi_2"*" [deg]", 
                framestyle=:box, labels=false,
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                markerstrokewidth = 0,
                title="t = "*string(time_snapshot_in_Myr)*" Myr")

    # scatter!(plt, [0], [0], label=false, color=:red)

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/snapshot_angles_stream_t_"*string(time_snapshot_in_Myr)*"_Myr.pdf"
    savefig(plt, namefile_pdf)


    s = 2.0

    plt = scatter(tab_Lz_wrt_host_in_kpckms_unbound[1:Nub] ./ 10^2, tab_E_wrt_host_in_kms_unbound[1:Nub] ./ 10^4,
                xlabel=L"L_z \ [10^2 \ \mathrm{kpc}\ \mathrm{km}/\mathrm{s}]", ylabel=L"E \ [10^4 \ (\mathrm{km}/\mathrm{s})^2]", 
                framestyle=:box, labels=false, 
                markerstrokewidth = 0,
                # aspect_ratio=1, 
                size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                # xticks=-3:0.5:3, yticks=-3:0.5:3,
                title="t = "*string(time_snapshot_in_Myr)*" Myr")

    display(plt)
    readline()

    mkpath(path_to_plot * "fig/"*name_model*"/")
    namefile_pdf = path_to_plot * "fig/"*name_model*"/snapshot_IOM_t_"*string(time_snapshot_in_Myr)*"_Myr.pdf"
    savefig(plt, namefile_pdf)

end

treat_data()