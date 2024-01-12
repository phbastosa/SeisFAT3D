#ifndef JARVIS_CONST_H
#define JARVIS_CONST_H
#include "jarvis_host_math_0_meat.hpp"
namespace jarvis
{
    template <typename eT>
    class jarvisConst
    {
    public:
        static const eT pi;          //*< ratio of any circle's circumference to its diameter
        static const eT e;           //*< base of the natural logarithm
        static const eT euler;       //*< Euler's constant, aka Euler-Mascheroni constant
        static const eT gratio;      //*< golden ratio
        static const eT sqrt2;       //*< square root of 2
        static const eT sqrt2pi;     //*< square root of 2*pi
        static const eT log_sqrt2pi; //*< log of square root of 2*pi
        static const eT eps;         //*< the difference between 1 and the least value greater than 1 that is representable
        //
        static const eT m_u;       //*< atomic mass constant (in kg)
        static const eT N_A;       //*< Avogadro constant
        static const eT k;         //*< Boltzmann constant (in joules per kelvin)
        static const eT k_evk;     //*< Boltzmann constant (in eV/K)
        static const eT a_0;       //*< Bohr radius (in meters)
        static const eT mu_B;      //*< Bohr magneton
        static const eT Z_0;       //*< characteristic impedance of vacuum (in ohms)
        static const eT G_0;       //*< conductance quantum (in siemens)
        static const eT k_e;       //*< Coulomb's constant (in meters per farad)
        static const eT eps_0;     //*< electric constant (in farads per meter)
        static const eT m_e;       //*< electron mass (in kg)
        static const eT eV;        //*< electron volt (in joules)
        static const eT ec;        //*< elementary charge (in coulombs)
        static const eT F;         //*< Faraday constant (in coulombs)
        static const eT alpha;     //*< fine-structure constant
        static const eT alpha_inv; //*< inverse fine-structure constant
        static const eT K_J;       //*< Josephson constant
        static const eT mu_0;      //*< magnetic constant (in henries per meter)
        static const eT phi_0;     //*< magnetic flux quantum (in webers)
        static const eT R;         //*< molar gas constant (in joules per mole kelvin)
        static const eT G;         //*< Newtonian constant of gravitation (in newton square meters per kilogram squared)
        static const eT h;         //*< Planck constant (in joule seconds)
        static const eT h_bar;     //*< Planck constant over 2 pi, aka reduced Planck constant (in joule seconds)
        static const eT m_p;       //*< proton mass (in kg)
        static const eT R_inf;     //*< Rydberg constant (in reciprocal meters)
        static const eT c_0;       //*< speed of light in vacuum (in meters per second)
        static const eT sigma;     //*< Stefan-Boltzmann constant
        static const eT R_k;       //*< von Klitzing constant (in ohms)
        static const eT b;         //*< Wien wavelength displacement law constant
        static const eT eps0;      //*< Vacuum dielectric constant(in F/m)

        // FOR CUDA
    public:
        static const int block_size;
    };

    // the long lengths of the constants are for future support of "long double"
    // and any smart compiler that does high-precision computation at compile-time

    template <typename eT>
    const eT jarvisConst<eT>::pi = eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);
    template <typename eT>
    const eT jarvisConst<eT>::e = eT(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274);
    template <typename eT>
    const eT jarvisConst<eT>::euler = eT(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495);
    template <typename eT>
    const eT jarvisConst<eT>::gratio = eT(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374);
    template <typename eT>
    const eT jarvisConst<eT>::sqrt2 = eT(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727);
    template <typename eT>
    const eT jarvisConst<eT>::sqrt2pi = eT(2.5066282746310005024157652848110452530069867406099383166299235763422936546078419749465958383780572661);
    template <typename eT>
    const eT jarvisConst<eT>::log_sqrt2pi = eT(0.9189385332046727417803297364056176398613974736377834128171515404827656959272603976947432986359541976);
    template <typename eT>
    const eT jarvisConst<eT>::m_u = eT(1.66053906660e-27);
    template <typename eT>
    const eT jarvisConst<eT>::N_A = eT(6.02214076e23);
    template <typename eT>
    const eT jarvisConst<eT>::k = eT(1.380649e-23);
    template <typename eT>
    const eT jarvisConst<eT>::k_evk = eT(8.617333262e-5);
    template <typename eT>
    const eT jarvisConst<eT>::a_0 = eT(5.29177210903e-11);
    template <typename eT>
    const eT jarvisConst<eT>::mu_B = eT(9.2740100783e-24);
    template <typename eT>
    const eT jarvisConst<eT>::Z_0 = eT(376.730313668);
    template <typename eT>
    const eT jarvisConst<eT>::G_0 = eT(7.748091729e-5);
    template <typename eT>
    const eT jarvisConst<eT>::k_e = eT(8.9875517923e9);
    template <typename eT>
    const eT jarvisConst<eT>::eps_0 = eT(8.8541878128e-12);
    template <typename eT>
    const eT jarvisConst<eT>::m_e = eT(9.1093837015e-31);
    template <typename eT>
    const eT jarvisConst<eT>::eV = eT(1.602176634e-19);
    template <typename eT>
    const eT jarvisConst<eT>::ec = eT(1.602176634e-19);
    template <typename eT>
    const eT jarvisConst<eT>::F = eT(96485.33212);
    template <typename eT>
    const eT jarvisConst<eT>::alpha = eT(7.2973525693e-3);
    template <typename eT>
    const eT jarvisConst<eT>::alpha_inv = eT(137.035999084);
    template <typename eT>
    const eT jarvisConst<eT>::K_J = eT(483597.8484e9);
    template <typename eT>
    const eT jarvisConst<eT>::mu_0 = eT(1.25663706212e-6);
    template <typename eT>
    const eT jarvisConst<eT>::phi_0 = eT(2.067833848e-15);
    template <typename eT>
    const eT jarvisConst<eT>::R = eT(8.314462618);
    template <typename eT>
    const eT jarvisConst<eT>::G = eT(6.67430e-11);
    template <typename eT>
    const eT jarvisConst<eT>::h = eT(6.62607015e-34);
    template <typename eT>
    const eT jarvisConst<eT>::h_bar = eT(1.054571817e-34);
    template <typename eT>
    const eT jarvisConst<eT>::m_p = eT(1.67262192369e-27);
    template <typename eT>
    const eT jarvisConst<eT>::R_inf = eT(10973731.568160);
    template <typename eT>
    const eT jarvisConst<eT>::c_0 = eT(299792458.0);
    template <typename eT>
    const eT jarvisConst<eT>::sigma = eT(5.670374419e-8);
    template <typename eT>
    const eT jarvisConst<eT>::R_k = eT(25812.80745);
    template <typename eT>
    const eT jarvisConst<eT>::b = eT(2.897771955e-3);
    template <typename eT>
    const eT jarvisConst<eT>::eps0 = eT(8.854187817e-12);
    // FOR CUDA
    template <typename eT>
    const int jarvisConst<eT>::block_size = 128; // 256 best for CUDA V100

    typedef jarvisConst<float> fjarvis_const;
    typedef jarvisConst<double> jarvis_const;
}
#endif