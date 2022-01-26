#ifndef MATH_INTEGRAL_H
#define MATH_INTEGRAL_H

// mohan add 2021-04-03
namespace ModuleBase
{

class Integral
{

	public:

    Integral();
    ~Integral();

    /**
     * @brief simpson integral.
     * 
     * @param mesh [in] number of grid points (should be odd)
     * @param func [in] function to be integrated
     * @param rab [in] a list of interval
     * @param asum [out] the final integral value
     */
    static void Simpson_Integral  // Peize Lin accelerate 2017-10-02
    (    
        const int mesh,
        const double * const func,
        const double * const rab,
        double &asum
    );

    /**
     * @brief simpson integral. 
     * 
     * @param mesh [in] number of grid points (should be odd)
     * @param func [in] function to be integrated
     * @param dr [in] interval
     * @param asum [out] the final integral value
     */
	static void Simpson_Integral  //Peize Lin accelerate 2017-10-02
	(
		const int mesh,
		const double * const func,
		const double dr,
		double &asum
	);

    /**
     * @brief simpson integral. 
     *  
     * 
     * @param mesh [in] number of grid points 
     * @param func [in] function to be integrated
     * @param rab [in] a list of interval
     * @param asum [out] a list of integral value. asum[i] = integral from 0 to i. The max index of asum is an even (mesh-1 or mesh).
     */
    static void Simpson_Integral_0toall  //Peize Lin add 2016-02-14.
    (
        const int mesh,
        const double * const func,
        const double * const rab,
        double * const asum
    );

    /**
     * @brief simpson integral. 
     * 
     * @param mesh [in] number of grid points 
     * @param func [in] function to be integrated 
     * @param rab [in] r(i) * dr(i)/di * di or (b-a)/2n 
     * @param asum [out] a list of integral value. sum[i] = integral from i to mesh-1. 
     */
    static void Simpson_Integral_alltoinf //Peize Lin add 2016-02-14.
    (
        const int mesh,
        const double * const func,
        const double * const rab,
        double * const asum
    );     

};

}
#endif
