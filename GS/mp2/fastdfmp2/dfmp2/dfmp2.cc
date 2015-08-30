#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libthce/lreri.h>
#include <libthce/thcew.h>
#include <libthce/laplace.h>
#include <libthce/thce.h>
#include <libciomr/libciomr.h>
#include <libplugin/plugin.h>
#include <libqt/qt.h>
#include <boost/shared_ptr.hpp>
#include <libmints/typedefs.h>
#include <lib3index/dftensor.h>
#include <lib3index/cholesky.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

INIT_PLUGIN

namespace psi {

extern "C" int read_options(std::string name, Options &options)
{
    if (name == "DFMP2"|| options.read_globals()) {

        /*- Print level -*/
        options.add_int("PRINT", 1);
        /*- Debug level -*/
        options.add_int("DEBUG", 0);
        /*- Bench level -*/
        options.add_int("BENCH", 0);
        options.add_double("CHOLESKY_TOLERANCE", 1e-14);

    }

    return true;
}

extern "C" PsiReturnType dfmp2(Options &options)
{
    tstart();

    outfile->Printf( "  ==> DF-MP2 <==\n\n");

    // => Setup and Sizing <= //

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<BasisSet> primary = wfn->basisset();

    boost::shared_ptr<BasisSetParser> parser (new Gaussian94BasisSetParser());
boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_orbital(primary->molecule(), "DF_BASIS_MP2",options.get_str("DF_BASIS_MP2"));
    
    boost::shared_ptr<Matrix> Caocc = wfn->Ca_subset("AO", "ACTIVE_OCC"); 
    boost::shared_ptr<Matrix> Cavir = wfn->Ca_subset("AO", "ACTIVE_VIR"); 
    boost::shared_ptr<Matrix> Cpq   = wfn->Ca_subset("AO", "ACTIVE");
    
    boost::shared_ptr<Vector> eps_aocc = wfn->epsilon_a_subset("AO", "ACTIVE_OCC"); 
    boost::shared_ptr<Vector> eps_avir = wfn->epsilon_a_subset("AO", "ACTIVE_VIR"); 

    int no = eps_aocc->dimpi()[0];   
    int nv = eps_avir->dimpi()[0];   
    int nn = primary->nbf();
    int nQ = auxiliary->nbf();

    long int memory = Process::environment.get_memory();
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    outfile->Printf("   => Sizing <=\n\n");
    outfile->Printf("    Memory  = %11zu MB\n", memory / (1024L * 1024L));
    outfile->Printf("    Threads = %11d\n", nthreads);
    outfile->Printf("    no      = %11d\n", no);
    outfile->Printf("    nv      = %11d\n", nv);
    outfile->Printf("    nn      = %11d\n", nn);
    outfile->Printf("    nQ      = %11d\n", nQ);
    outfile->Printf("\n");
    outfile->Printf("   => Molecule <=\n\n");
    primary->molecule()->print();
    outfile->Printf("   => Primary Basis <=\n\n");
    primary->print();
    outfile->Printf("   => Auxiliary Basis <=\n\n");
 
    // => Setup DF Integrals <= //

    //The constructor for building density fitted ints
    int nfvir = wfn->frzvpi().sum();
    int nfocc = wfn->frzcpi().sum();
    int naocc = wfn->nalphapi().sum() - nfocc;
    int navir = wfn->nmopi().sum() - naocc - nfocc - nfvir;
    int nmo   = naocc+navir;
    outfile->Printf("\nOrbital Information\n");
    outfile->Printf("\n nfocc = %d  nfvir = %d\n", nfocc, nfvir);
    outfile->Printf("\n naocc = %d  navir = %d\n", naocc, navir);
    outfile->Printf("\n NMO = %d", nmo);
    
    //This is active in the vain of SR methods.  No frozen orbitals.  
    //Contains Bpq where pq are any arbitrary integral
    //If I want to include all orbitals, just pass "ALL" and make sure I do not subtract frozen
    boost::shared_ptr<DFERI> df = DFERI::build(primary,auxiliary,options, wfn);
    df->add_pair_space("B", "ACTIVE_ALL", "ACTIVE_ALL");

    df->set_memory(memory / 8L);
    df->print_header();
    df->compute();

    boost::shared_ptr<Tensor> B = df->ints()["B"];
    df.reset();
    

    //FILE* Bf = B->file_pointer();


    // => DFMP2 Energy Evaluation <= //

    // => Setup disk blocks and integral buffers <= //

    long int doubles = memory / 8L;
    long int nvQ = nv * (long int) nQ;
    
    long int max_o = doubles / (2L * nvQ);
    max_o = (max_o > no ? no : max_o);
   
    //Different ways to store matrix - first is for DGEMM later and second is for storing matrix
    boost::shared_ptr<Matrix> pqB(new Matrix("pqB",(naocc+navir)*(naocc+navir),nQ )); 
    boost::shared_ptr<Matrix> Bpq(new Matrix("Bpq",(naocc+navir),(naocc+navir)*nQ )); 

    double** Bpqp = Bpq->pointer();
    double** pqBp = pqB->pointer();
    FILE* Bf = B->file_pointer();
    
  
    //for(int p = 0; p < naocc+navir; p++){
    
    //This assumes that the code can all fit in core.  
    //Rob's original code does this in blocks.  I could use that if I wanted.  
    //The code is arranged in the file in order of fastest index: naux, nmo, nmo
    fseek(Bf,0, SEEK_SET);
    fread(Bpqp[0], sizeof(double),nQ*(naocc+navir)*(naocc+navir), Bf); 
    //EMP2_corr = \sum_{ijab} (ia|jb) * (2(ia|jb) - (ib|ja))/(occ - vir)
    
    SharedMatrix pqrs(new Matrix("pqrs", (naocc+navir)*(naocc+navir), (naocc+navir)*(naocc+navir)));
    SharedMatrix pqrs_gemm(new Matrix("pqrs_gemm", (naocc+navir)*(naocc+navir), (naocc+navir)*(naocc+navir)));
    //This was used as a test to see if I read the integrals in correctly.  I did!
    //See if I can generate CD integrals

    outfile->Printf("\n Generating CD Integrals");
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary, primary, primary, primary));
    //This is the tolerance for CD (very low = exact)
    double tol_cd = options.get_double("CHOLESKY_TOLERANCE");
    boost::shared_ptr<CholeskyERI> Ch (new CholeskyERI(boost::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0 ,tol_cd, Process::environment.get_memory()));
    //Performs a Cholesky decomposition.  
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix aoeri(new Matrix("aoeri", nn*nn, nn*nn));
    aoeri = mints->ao_eri();
    //aoeri->print_to_mathematica();

    Ch->choleskify();
    int nL = Ch->Q();
    outfile->Printf("\n Number of cholesky vectors %d to satisfy %20.12f tolerance\n", nL,tol_cd);
    //SharedMatrix D = wfn->Da();
    
    //This is not necessary for MP2, but this will perform a cholesky decomp of any matrix that you give thats semi-positive definite
    //boost::shared_ptr<CholeskyLocal> ChM(new CholeskyLocal(Caocc, 1e-12, Process::environment.get_memory()));
    ////This was to see if cholesky can decompose AO Density matrix
    //ChM->choleskify();
    //SharedMatrix Clocalocc = ChM->L();
    //boost::shared_ptr<CholeskyLocal> ChMvir(new CholeskyLocal(Cavir, 1e-12, Process::environment.get_memory()));
    //ChMvir->choleskify();
    //SharedMatrix Clocalvir = ChMvir->L();

    //SharedMatrix Cfull(new Matrix("Cfull", nmo, nmo));
    //for(int i = 0; i < naocc; i++){
    //   for(int j = 0; j < nn; j++){
    //      Cfull->set(i,j,Clocalocc->get(i,j));
    //   }
    //}
    //for(int a = 0; a < navir; a++){
    //   for(int j = 0; j < nn; j++){
    //      Cfull->set(a+naocc,j,Clocalvir->get(a,j));  
    //   }
    //}
    //Cpq->copy(Cfull);
    //Cpq->print_to_mathematica();

    
    
    
    
    SharedMatrix Lao = Ch->L();
    //Lao->print_to_mathematica();
    SharedMatrix L(new Matrix("Lmo", nL, (naocc+navir)*(naocc+navir)));
    //This step transforms the Lao matrix to the Lmo.  should optimize.  use GEMM
    //Use a density matrix to transform both indices 
    for(int l = 0; l < nL; l++){
       for(int p = 0; p < naocc+navir; p++){
          for(int q = 0; q < naocc+navir; q++){
             for(int mu = 0; mu < nn; mu++){
                for(int nu = 0; nu < nn; nu++){
                   L->add(l,p*(naocc+navir)+q,Lao->get(l,mu*(nn)+nu)*Cpq->get(mu,p)*Cpq->get(nu,q)); 
                }
             }
          }
       }
    }
    //L->print_to_mathematica();
    
    double val = 0.0;
    double valk = 0.0;
    outfile->Printf("\nI am debuggin =:(");
    outfile->Printf("\n Printing index followed by value");
    for(int p = 0; p < (naocc+navir); p++){
      for(int q = 0; q < (naocc+navir); q++){
         for(int r = 0; r < (naocc+navir); r++){
            for(int s = 0; s < (naocc+navir); s++){
                for(int B = 0; B < nL; B++){
                    int qB = q*nQ + B;
                    int sB = s*nQ + B;
                    //val+=Bpq->get(p,qB)*Bpq->get(r,sB);
                    //This is done to resort my integrals into pq by B.  
                    //Makes GEMM call very trival(contract the B index)
                    pqB->set(p*(naocc+navir)+q,B,Bpq->get(p,qB));
    //                val+=L->get(B,p*(naocc+navir)+q)*L->get(B,r*(naocc+navir) + s); 
                    
                }
                
                outfile->Printf("\n p = %d q = %d r = %d s = %d %20.12f", p, q, r, s, val);
                int pq = p*(naocc+navir) + q;
                int rs = r*(naocc+navir) + s;
                //pqrs_gemm->set(pq,rs,val);
     
                val = 0.0;
             }
          }
       }
    }
    pqrs_gemm->gemm('N','T',(naocc+navir)*(naocc+navir),(naocc+navir)*(naocc+navir),nQ,1.0,pqB,nQ,pqB,nQ,0.0,(naocc+navir)*(naocc+navir),0,0,0);
    pqrs->gemm('T','N',(naocc+navir)*(naocc+navir),(naocc+navir)*(naocc+navir),nL,1.0,L,(naocc+navir)*(naocc+navir),L,(naocc+navir)*(naocc+navir),0.0,(naocc+navir)*(naocc+navir),0,0,0);
    //pqrs->gemm('T','N',1.0,L,L,0.0);
    
    double mp2_energy_L = 0.0;
    double mp2_energy_D = 0.0;
    for(int i = 0; i < naocc; i++){
       for(int a = naocc; a < naocc+navir; a++){
          for(int j = 0; j < naocc; j++){
             for(int b = naocc; b < naocc+navir; b++){
               double denom = eps_aocc->get(i) + eps_aocc->get(j) - eps_avir->get(a-naocc) - eps_avir->get(b-naocc);
               //Calculate MP2_energy for both Cholesky and Density fitting
               mp2_energy_L+=(pqrs->get(i*(naocc+navir)+a,j*(naocc+navir) + b))*(2.0*pqrs->get(i*(naocc+navir)+a,j*(naocc+navir) + b) - pqrs->get(i*(naocc+navir) + b,j*(naocc+navir) + a))/(denom);
               mp2_energy_D+=(pqrs_gemm->get(i*(naocc+navir)+a,j*(naocc+navir) + b))*(2.0*pqrs_gemm->get(i*(naocc+navir)+a,j*(naocc+navir) + b) - pqrs_gemm->get(i*(naocc+navir) + b,j*(naocc+navir) + a))/(denom);
             }
          }
       }
    }
    //pqrs->print_to_mathematica();    
    boost::shared_ptr<DFTensor> df3 (new DFTensor(primary, auxiliary, Cpq, naocc, navir, naocc, navir, options));
    SharedMatrix aopqrs(new Matrix("aopqrs", nn*nn, nn*nn));
    aopqrs = df3->Imo();
    //aopqrs->print_to_mathematica();

    outfile->Printf("\nMP2_energy = %20.12f for laplace", mp2_energy_L);
    outfile->Printf("\nMP2_energy = %20.12f for density fitting", mp2_energy_D);
    outfile->Printf("\nERROR = %20.12f", fabs(mp2_energy_D - mp2_energy_L));
  
    double* eop = eps_aocc->pointer();
    double* evp = eps_avir->pointer();

    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    boost::shared_ptr<LaplaceDenom> laplace(new LaplaceDenom::LaplaceDenom(eps_aocc, eps_avir, 1e-6, 0.0, 2));
    laplace->compute("pi_i", "pi_a");
    //thce->add_tensor("pi_i", laplace->tau_occ());
    //thce->add_tensor("pi_a", laplace->tau_vir());
    boost::shared_ptr<Tensor> tau_occ = laplace->tau_occ();
    boost::shared_ptr<Tensor> tau_vir = laplace->tau_vir();
    tau_occ->print();
    //(*thce)["pi_i"]->dimensions()[0] = "nw";
    //(*thce)["pi_i"]->dimensions()[1] = "naocc";
    //(*thce)["pi_a"]->dimensions()[0] = "nw";
    //(*thce)["pi_a"]->dimensions()[1] = "navir";
    //thce->new_dimension("nw", laplace->npoints());
    tau_vir->print();
    double* tau_occp = tau_occ->pointer();
    double* tau_virp = tau_vir->pointer();
    double denom_test = 0.0;
    double Delta = 0.0;
    for(int i = 0; i < no; i++){
       for(int j = 0; j < no; j++){
          for(int a = 0; a < nv; a++){
             for(int b = 0; b < nv; b++){
                Delta = eps_aocc->get(i) + eps_aocc->get(j) - eps_avir->get(a) - eps_avir->get(b);
                denom_test += (1.0 - exp(-1*Delta*Delta))/(Delta);
             }
          }
       }
    }
    outfile->Printf("Denom_test = %20.12f", denom_test);
    denom_test = 0.0;
    double PI = 0.0;
    double PA = 0.0;
    double PJ = 0.0;
    double PB = 0.0;
    for(int i = 0; i < no; i++){
       for(int omega = 0; omega < tau_occ->active_sizes()[0]; omega++){
             int om = tau_occ->active_sizes()[1];
             outfile->Printf("\ni =  %d  tau_occ = %20.12f   PI = %20.12f\n",i, *(tau_occp + om*omega + i), PI);
             PI+=*(tau_occp + om*omega+ i);
       } 
       for(int j = 0; j < no; j++){
          for(int omega = 0; omega < tau_occ->active_sizes()[0]; omega++){
             int om = tau_occ->active_sizes()[1];
             PJ+=*(tau_occp + om*omega+ j);
          } 
          for(int a = 0; a < nv; a++){
             for(int omega = 0; omega < tau_vir->active_sizes()[0]; omega++){
                int om = tau_vir->active_sizes()[1];  
                PA+=*(tau_virp + om*omega + a);
             }
             for(int b = 0; b < nv; b++){
                for(int omega = 0; omega < tau_vir->active_sizes()[0]; omega++){
                   int om = tau_vir->active_sizes()[1];  
                   PB+=*(tau_virp + om*omega + b);
             }
                Delta = eps_aocc->get(i) + eps_aocc->get(j) - eps_avir->get(a) - eps_avir->get(b);
                denom_test += (1.0 - exp(-1*Delta*Delta))*PI*PJ*PA*PB;
             PI =0.0;
             PJ =0.0;
             PA =0.0;
             PB =0.0;
             }
          }
       }
    }
    outfile->Printf("\nDenom_test = %20.12f 'n", denom_test);
    // => Form energy contributions <= //

    //for (int istart = 0; istart < no; istart += max_o) {
    //    int ni = (istart + max_o >= no ? no - istart : max_o);


    //    for (int jstart = 0; jstart < no; jstart += max_o) {
    //        int nj = (jstart + max_o >= no ? no - jstart : max_o);

    //        if (jstart > istart) break;

    //        if (istart == jstart) {
    //            Bjb->copy(Bia);
    //        } else {
    //            fseek(Bf,jstart * nvQ * sizeof(double), SEEK_SET);
    //            fread(Bjbp[0], sizeof(double), nj * nvQ, Bf); 
    //        }
    //        
    //        #pragma omp parallel for reduction(+: E_MP2J, E_MP2K)
    //        for (int ijrel = 0; ijrel < ni *nj; ijrel++) {        

    //            int irel = ijrel / nj;
    //            int jrel = ijrel % nj;

    //            int i = irel + istart; 
    //            int j = jrel + jstart; 

    //            if (j > i) continue;

    //            int thread = 0;
    //            #ifdef _OPENMP
    //                thread = omp_get_thread_num();
    //            #endif

    //            double** Iabp = Iab[thread]->pointer();

    //            double perm = (i == j ? 1.0 : 2.0);
    //        
    //            C_DGEMM('N','T',nv,nv,nQ,1.0,Biap[irel],nQ,Bjbp[jrel],nQ,0.0,Iabp[0],nv);

    //            for (int a = 0; a < nv; a++) {
    //                for (int b = 0; b < nv; b++) {
    //                    double iajb = Iabp[a][b];
    //                    double ibja = Iabp[b][a];
    //                    double D = perm / (eop[i] + eop[j] - evp[a] - evp[b]);
    //                    E_MP2J += 2.0 * iajb * iajb * D;
    //                    E_MP2K -= 1.0 * iajb * ibja * D;
    //                }
    //            }
    //        }
    //    }
    //}
    //    

    //double E_MP2 = E_MP2J + E_MP2K;

    //outfile->Printf( "    @DF-MP2 Correlation Energy: %24.16f\n", E_MP2);
   //*/
    tstop();

    return Success;
}

} // end namespaces
