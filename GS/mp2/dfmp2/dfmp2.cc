/*
 *@BEGIN LICENSE
 *
 * dfmp2 by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/wavefunction.h>
#include <libmints/matrix.h>
#include <libpsio/psio.h>
#include <boost/multi_array.hpp>
#include <lib3index/3index.h>
#include <libqt/qt.h>
#include <libthce/lreri.h>
#include <libthce/thce.h>

INIT_PLUGIN

namespace psi{ namespace dfmp2 {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "DFMP2"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
    }

    return true;
}

extern "C"
PsiReturnType dfmp2(Options &options)
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");
    boost::shared_ptr<BasisSet> ribasis = BasisSet::construct(parser, molecule, "DF_BASIS_MP2");
    int naux = ribasis->nbf();
    int nbf[] ={ aoBasis->nbf()  };
    int nbfA = aoBasis->nbf();
    
    outfile->Printf("\nnaux = %d  nbf = %d", naux, nbfA);
    SharedMatrix Caocc = wfn->Ca_subset("AO", "ACTIVE_OCC");
    SharedMatrix Cavir = wfn->Ca_subset("AO", "ACTIVE_VIR");

    int naocc = Caocc->colspi()[0];
    int navir = Cavir->colspi()[0];
    int maxQ  = ribasis->max_function_per_shell();
   
    typedef boost::multi_array<double, 3> third;
    typedef boost::multi_array<double, 4> fourth;
    third Amn(boost::extents[naux][nbfA][nbfA]);
    third Ami(boost::extents[naux][nbfA][naocc]);
    third Aia(boost::extents[naux][naocc][navir]);
    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (ribasis, BasisSet::zero_ao_basis_set(), aoBasis, aoBasis));
    boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
    const double *buffer = eri->buffer();
    
    SharedMatrix AmnM(new Matrix("Amn", naux, nbfA*nbfA)); 
    SharedMatrix AmiM(new Matrix("Ami", naux, nbfA*naocc)); 
    SharedMatrix AiaM(new Matrix("Aia", naux, navir*naocc)); 
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    double nucrep = molecule->nuclear_repulsion_energy();
    outfile->Printf( "\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);

    // The matrix factory can create matrices of the correct dimensions...

    // Form the one-electron integral objects from the integral factory
    // Form the one-electron integral matrices from the matrix factory
    double** AmnMp = AmnM->pointer();
    for(int Qshell = 0; Qshell < ribasis->nshell(); Qshell++){
       for(int Mshell = 0; Mshell < aoBasis->nshell(); Mshell++){
          for(int Nshell = 0; Nshell < aoBasis->nshell(); Nshell++){
             eri->compute_shell(Qshell, 0, Mshell, Nshell);
             int nQ = ribasis->shell(Qshell).nfunction();
             int nM = aoBasis->shell(Mshell).nfunction();
             int nN = aoBasis->shell(Nshell).nfunction();
             int oQ = ribasis->shell(Qshell).function_index();
             int oM = aoBasis->shell(Mshell).function_index();
             int oN = aoBasis->shell(Nshell).function_index();
             
             for(int Q = 0, index = 0; Q < nQ; Q++){
                for(int M = 0; M < nM; M++){
                   for(int N = 0; N < nN; N++,index++){
                      AmnMp[Q+oQ][(M + oM)*nbfA + N + oN] = buffer[index];
                      AmnMp[Q+oQ][(N + oN)*nbfA + M + oM] = buffer[index];
                      //outfile->Printf("\nAmn[%d][%d][%d] = %8.6f", oQ,oM,oN,buffer[index]);
                      //outfile->Printf("\nAmn[%d][%d][%d]",Q + oQ,M + oM, N + oN );
                   
                   } 
                }
              }
          }
       }
    }  
    AmiM->gemm('N','N', naux*nbfA, naocc, nbfA, 1.0, AmnM, nbfA, Caocc, naocc, 0.0, naocc);

    
    for(int A = 0; A < naux; A++){
       C_DGEMM('T','N',naocc, navir, nbfA, 1.0, AmiM->pointer()[A],naocc, Cavir->pointer()[0], navir, 0.0, AiaM->pointer()[A], navir);
    }
    
    boost::shared_ptr<FittingMetric> metric(new FittingMetric(ribasis, true));
    metric->form_eig_inverse(1.0E-10);
    SharedMatrix Jm12 = metric->get_metric();
    third Qia(boost::extents[naux][naocc][navir]);
    SharedMatrix QiaM(new Matrix("Qia", naux, naocc*navir));

    for(int Q = 0; Q < naux; Q++){
       for(int i = 0; i < naocc; i++){
          for(int a = 0; a < navir; a++){
             Aia[Q][i][a] = AiaM->get(Q,i*navir + a); 
          }
       }
    }
    for(int Q = 0; Q < naux; Q++){
       for(int A = 0; A < naux; A++){
          for(int i = 0; i < naocc; i++){
             for(int a = 0; a < navir; a++){
                QiaM->add(Q,i*navir + a, Jm12->get(Q,A)*AiaM->get(A,i*navir+a));;
                Qia[Q][i][a] += Jm12->get(Q,A)*Aia[A][i][a];
             }
          }
       }
    }
    //outfile->Printf("Qia_sm[1][1][1] = %20.12f \ n Qia_bm[1][1][1] = %20.12f", QiaM->get(1,1*navir + 1), Qia[1][1][1]);
    //QiaM->zero();
    //QiaM->gemm('T','N', naux, naocc*navir, naux, 1.0, AiaM,naux, Jm12, naux,0.0,navir);
     
    SharedVector eps_avir = wfn->epsilon_a_subset("AO", "ACTIVE_VIR"); 
    SharedVector eps_aocc = wfn->epsilon_a_subset("AO", "ACTIVE_OCC"); 
    double MP2_energy = 0.0;
    SharedMatrix dfmo(new Matrix("DFMO",naocc*navir,naocc*navir));
    for(int i = 0; i < naocc; i++){
       for(int j = 0; j < naocc; j++){
              
              //C_DGEMM('T','N',navir, navir, naux, 1.0,QiaM->pointer()[i*navir],navir,QiaM->pointer()[j*navir],naux,0.0,Iab->pointer()[0],navir);
              for(int a = 0; a < navir; a++){
                 for(int b =0; b <navir; b++){
                    for(int Q = 0; Q < naux; Q++){
                         dfmo->add(i*navir + a,j*navir + b,QiaM->get(Q,i*navir + a)*QiaM->get(Q,j*navir + b));
                    }
                 }
              }
       }
    }
    for(int i = 0; i < naocc; i++){
       for(int j = 0; j < naocc; j++){
          for(int a = 0; a < navir; a++){
             for(int b = 0; b < navir; b++){
                       double iajb = dfmo->get(i*navir + a, j*navir + b); 
                       double ibja = dfmo->get(i*navir + b, j*navir + a); 
                       double denom = -1.0 / (eps_avir->get(a) + eps_avir->get(b) - eps_aocc->get(i) - eps_aocc->get(j));
                       MP2_energy+= 2.0*(iajb*iajb)*denom - iajb*ibja*denom;
             }
          }
       }
    }
    outfile->Printf("\n\t MP2_energy = %20.12f", MP2_energy);

             
          
    
    SharedMatrix C = wfn->Ca();
    boost::shared_ptr<DFTensor> DF(new DFTensor(aoBasis, ribasis, C, naocc, navir, naocc, navir, options));
    SharedMatrix Imo = DF->Idfmo();
  


      /*    for(int ijrel = 0; ijrel < naocc*naocc; ijrel++){
             
             int irel = ijrel / naocc;
             int jrel = ijrel % naocc;
 
             int i = irel + i;
             int j = jrel + j;

             if(j > i) continue;
  
             double perm = (i ==j ? 1.0 : 2.0);
             C_DGEMM('N','T', navir, navir, naux, 1.0, Bia->pointer()[irel], naux, Bjb->pointer()[jrel], naux, 0.0,Iab->pointer()[0], navir);

             for(int a = 0; a < navir; a++){
                for(int b = 0; b < navir; b++){
                    double iajb = Iab->get(a,b);
                    double ibja = Iab->get(b,a);

                   double D = perm / (eps_aocc->get(i) + eps_aocc->get(j) - eps_avir->get(a) - eps_avir->get(b));
                   robplmp2+=2.0 * iajb * iajb *D;
                   robplmp2-=1.0 * iajb*ibja*D;
             
                }
            }
         }*/
      //}
    //} 
    
   /* SharedMatrix Imo= DF->Idfmo();
    //SharedMatrix Imo= DF->Imo();
    Imo->print();
    MP2_energy = 0.0;
    for(int i = 0; i < naocc; i++){
       for(int j = 0; j < naocc; j++){
          for(int a = naocc; a < naocc + navir; a++){
             for(int b = naocc; b < naocc + navir; b++){
                double perm_factor = (i == j ? 1.0 : 2.0);
                double denom = perm_factor / (eps_avir->get(a - naocc) + eps_avir->get(b - naocc) - eps_aocc->get(i) - eps_aocc->get(j));
                int ia = i*navir + a;
                int jb = j*navir + b;
                int ib = i*navir + b;
                int ja = j*navir + a;
                double iajb = Imo->get(ia,jb);
                double ibja = Imo->get(ib,ja);
                MP2_energy+=iajb*iajb*denom;
                MP2_energy-=(iajb*ibja*denom);
             }
          }
       }
    }
    outfile->Printf("\n\n\t Energy with Imo: %20.12f", MP2_energy);
    */ 
   
    
    return Success;
}

}} // End Namespaces
