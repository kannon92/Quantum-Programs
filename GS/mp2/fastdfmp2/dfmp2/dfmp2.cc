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
    boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary->molecule(), "DF_BASIS_MP2");
    
    boost::shared_ptr<Matrix> Caocc = wfn->Ca_subset("AO", "ACTIVE_OCC"); 
    boost::shared_ptr<Matrix> Cavir = wfn->Ca_subset("AO", "ACTIVE_VIR"); 
    boost::shared_ptr<Matrix> Cpq   = wfn->Ca();
    
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
    auxiliary->print();
 
    // => Setup DF Integrals <= //

    boost::shared_ptr<DFERI> df = DFERI::build(primary,auxiliary,options, wfn);
    int nfvir = wfn->frzvpi().sum();
    int nfocc = wfn->frzcpi().sum();
    int naocc = wfn->nalphapi().sum() - nfocc;
    int navir = wfn->nmopi().sum() - naocc - nfocc - nfvir;
    
    //This is active in the vain of SR methods.  No frozen orbitals.  
    //Contains Bpq where pq are any arbitrary integral
    df->add_pair_space("B", "ACTIVE_ALL", "ACTIVE_ALL");

    df->set_memory(memory / 8L);
    df->print_header();
    df->compute();

    //boost::shared_ptr<Tensor> pqB =CoreTensor::build("B", "NAUX",nQ, "NALL", naocc+navir, "NALL", naocc+navir);
    boost::shared_ptr<Tensor> B = df->ints()["B"];
    df.reset();
    //thce->print();
    //pqB->print();
    //B->swap_in(true);
    

    //FILE* Bf = B->file_pointer();


    // => DFMP2 Energy Evaluation <= //

    // => Setup disk blocks and integral buffers <= //

    long int doubles = memory / 8L;
    long int nvQ = nv * (long int) nQ;
    
    long int max_o = doubles / (2L * nvQ);
    max_o = (max_o > no ? no : max_o);
   
    boost::shared_ptr<Matrix> Bpq(new Matrix("pqB",(naocc+navir),(naocc+navir)*nQ )); 
    boost::shared_ptr<Matrix> Brs(new Matrix("rsB",(naocc+navir),(naocc+navir)*nQ )); 

    double** Bpqp = Bpq->pointer();
    double** Brsp = Brs->pointer();
    FILE* Bf = B->file_pointer();
    B->print();
    
  
    //This dimension because data is layed out (mo*mo) by Q
    //This is the part of the code I am not confident in.  
    //How exactly do I read into from the disk.  
    //These are layed out in slowest to fastest by p by q by naux
    //for(int p = 0; p < naocc+navir; p++){
    fseek(Bf,(naocc+navir)*(naocc+navir)*nQ*sizeof(double), SEEK_SET);
    fread(Bpqp[0], sizeof(double),nQ*(naocc+navir)*(naocc+navir), Bf); 
    //}
    Bpq->print();
    //Bpq->print();
  
    //These are generated using lib3index.  This could be used as a way to test my code.  
    //One major problem with lib3index is not this code is not parallezied like Rob's tensor code.  
    
    boost::shared_ptr<DFTensor> DF(new DFTensor(primary, auxiliary, Cpq, naocc, navir, naocc, navir, options));
    //The Qpq term - MO int with fitted folded in
    SharedMatrix Qpq = DF->Qmo();
    Qpq->print();
    //Qpq * Qrs = pq by rs) - Fully two ints 
    //SharedMatrix dfmo = DF->Idfmo();
    //dfmo->print();
     

    /*std::vector<boost::shared_ptr<Matrix> > Iab;
    for (int t = 0; t < nthreads; t++) {
        Iab.push_back(boost::shared_ptr<Matrix>(new Matrix("Iab", nv, nv)));
    }

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

    for (int istart = 0; istart < no; istart += max_o) {
        int ni = (istart + max_o >= no ? no - istart : max_o);


        for (int jstart = 0; jstart < no; jstart += max_o) {
            int nj = (jstart + max_o >= no ? no - jstart : max_o);

            if (jstart > istart) break;

            if (istart == jstart) {
                Bjb->copy(Bia);
            } else {
                fseek(Bf,jstart * nvQ * sizeof(double), SEEK_SET);
                fread(Bjbp[0], sizeof(double), nj * nvQ, Bf); 
            }
            
            #pragma omp parallel for reduction(+: E_MP2J, E_MP2K)
            for (int ijrel = 0; ijrel < ni *nj; ijrel++) {        

                int irel = ijrel / nj;
                int jrel = ijrel % nj;

                int i = irel + istart; 
                int j = jrel + jstart; 

                if (j > i) continue;

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Iabp = Iab[thread]->pointer();

                double perm = (i == j ? 1.0 : 2.0);
            
                C_DGEMM('N','T',nv,nv,nQ,1.0,Biap[irel],nQ,Bjbp[jrel],nQ,0.0,Iabp[0],nv);

                for (int a = 0; a < nv; a++) {
                    for (int b = 0; b < nv; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double D = perm / (eop[i] + eop[j] - evp[a] - evp[b]);
                        E_MP2J += 2.0 * iajb * iajb * D;
                        E_MP2K -= 1.0 * iajb * ibja * D;
                    }
                }
            }
        }
    }
        

    double E_MP2 = E_MP2J + E_MP2K;

    outfile->Printf( "    @DF-MP2 Correlation Energy: %24.16f\n", E_MP2);
  e */
    tstop();

    return Success;
}

} // end namespaces