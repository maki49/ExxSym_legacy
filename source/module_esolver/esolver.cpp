#include "esolver.h"
#include "esolver_ks_pw.h"
#include "esolver_ks_lcao.h"
#include "esolver_of.h"

namespace ModuleESolver
{
void ESolver:: printag()
{
    std::cout<<tag<<std::endl;
}


//Some API to operate E_Solver
void init_esolver(ESolver *&p_esolver, const string use_esol)
{
     if(use_esol == "ksdft_pw")
     {
         p_esolver = new ESolver_KS_PW();
     }
    else if(use_esol == "ksdft_lcao")
     {
         p_esolver = new ESolver_KS_LCAO();
     }
     //  else if(use_esol == "sdft_pw")
    //  {
    //      p_esolver = new ESolver_KS_PW(true);
    //  }
    //  else if(use_esol == "ofdft")
    //  {
    //      p_esolver = new OFDFT();
    //  }    
}

void clean_esolver(ESolver *&pesolver)
{
    if(pesolver!=NULL) 
    {
        delete pesolver;
    }
}

}