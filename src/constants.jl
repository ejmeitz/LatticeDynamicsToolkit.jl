const periodic_table = PeriodicTable.elements

const lo_tol = 1e-5
const lo_sqtol = lo_tol^2
const lo_freqtol = lo_tol*1e-4

# Energy
const Hartree_to_eV = 27.21138602
const eV_to_Hartree = 1.0 / Hartree_to_eV
const eV_to_Joule=1.6021766208E-19
const Joule_to_eV=1.0/eV_to_Joule
const Hartree_to_Joule=Hartree_to_eV*eV_to_Joule
const Joule_to_Hartree=1.0/Hartree_to_Joule


# Planck Constant
const hbar_eV = 6.582119514E-16
const hbar_Hartree = hbar_eV*eV_to_Hartree
const hbar_Joule= hbar_eV* eV_to_Joule

# Length
const bohr_to_A = 0.52917721067
const A_to_bohr = 1.0 / bohr_to_A
const bohr_to_m=0.52917721067E-10
const m_to_bohr=1.0/bohr_to_m

# Mass (electorn mass unit / atomic mass unit)
const amu_to_kg = 1.660539040E-27
const emu_to_kg = 9.10938356E-31
const amu_to_emu = amu_to_kg/emu_to_kg
const emu_to_amu = 1.0/amu_to_emu

# Time
const time_au_to_s=hbar_Joule/Hartree_to_Joule
const time_s_to_au=1.0/time_au_to_s
const time_au_to_fs=1E15*time_au_to_s
const time_fs_to_au=1.0/time_au_to_fs

# Velocities
const velocity_au_to_ms=bohr_to_m/time_au_to_s
const velocity_ms_to_au=1.0/velocity_au_to_ms
const velocity_au_to_Afs=velocity_au_to_ms*1E-5
const velocity_Afs_to_au=1.0/velocity_au_to_Afs

const groupvel_Hartreebohr_to_ms=bohr_to_m/time_au_to_s
const groupvel_ms_to_Hartreebohr=1.0/groupvel_Hartreebohr_to_ms

const kB_eV = 8.6173303E-5
const kB_Hartree = kB_eV*eV_to_Hartree

const forceconstant_2nd_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^2)
const forceconstant_3rd_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^3)
const forceconstant_4th_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^4)

const forceconstant_2nd_HartreeBohr_to_eVA = 1.0 / forceconstant_2nd_eVA_to_HartreeBohr
const forceconstant_3rd_HartreeBohr_to_eVA = 1.0 / forceconstant_3rd_eVA_to_HartreeBohr
const forceconstant_4th_HartreeBohr_to_eVA = 1.0 / forceconstant_4th_eVA_to_HartreeBohr

const frequency_Hartree_to_THz=1e-12/(2*pi)/hbar_Hartree
const frequency_THz_to_Hartree=1.0/frequency_Hartree_to_THz
const frequency_Hartree_to_rads=1.0/hbar_Hartree


const _WELLDEFINED_SMALL_64 = [
    0.0,
    0.1, -0.1,
    0.1111111111111111, -0.1111111111111111,
    0.125, -0.125,
    0.1428571428571428, -0.1428571428571428,
    0.1666666666666667, -0.1666666666666667,
    0.2, -0.2,
    0.2222222222222222, -0.2222222222222222,
    0.25, -0.25,
    0.2857142857142857, -0.2857142857142857,
    0.3, -0.3,
    0.3333333333333333, -0.3333333333333333,
    0.375, -0.375,
    0.4, -0.4,
    0.4285714285714285, -0.4285714285714285,
    0.4444444444444444, -0.4444444444444444,
    0.5, -0.5,
    0.5555555555555556, -0.5555555555555556,
    0.5714285714285714, -0.5714285714285714,
    0.6, -0.6,
    0.625, -0.625,
    0.6666666666666667, -0.6666666666666667,
    0.7, -0.7,
    0.7071067811865475, -0.7071067811865475,   # 1/√2
    0.7142857142857143, -0.7142857142857143,
    0.75, -0.75,
    0.7777777777777778, -0.7777777777777778,
    0.8, -0.8,
    0.8333333333333334, -0.8333333333333334,
    0.8571428571428571, -0.8571428571428571,
    0.8660254037844386, -0.8660254037844386,   # √3/2
    0.875, -0.875,
    0.8888888888888888, -0.8888888888888888,
    0.9, -0.9,
    1.0, -1.0
]  # length = 69