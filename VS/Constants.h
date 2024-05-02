#ifndef CONSTANTS_H
#define CONSTANTS_H

class Constants
{
	public:
		static const double
			//----------- h2_mN = h_bar^2*c^2/amu (MeV fm^2) ----------------------------
			h2_mn,/*MeV fm**2*/
			//---------------------------------------------------------------------------
			// fundamental constants
			// (from Fundamental physics constants, UFN 173 (2003) p.339)
			hbar, // 10^(-22) MeV sec
			c,    // 10^(15) fm/sec
			mp,   // proton mass (MeV)
			mn,   // neutron mass (MeV)
			mN,   // (mp+mn)/2 (MeV)
			amu,  // atomic mass unit (MeV)
			h_mn, // hbar/mN (10^(23) fm^2/sec)
			e2,   // MeV*fm
			hc,   // hbar*c (MeV*fm)
			kB,   // Boltzmann constant, MeV/10^9K
			NA,   // Avagadro constant 10^23 mol^-1     
			zerro,
			small_add,
			small_e,
			large_val,
			degTorad,
			radTodeg,
			M_PI, 
			M_PI_2, 
			M_PI_4, 
			M_1_PI, 
			M_2_PI, 
			M_2_SQRTPI, M_SQRT2, M_SQRT1_2,
			M_1_3;
		static const bool 
			CondonShortley;
};
#endif