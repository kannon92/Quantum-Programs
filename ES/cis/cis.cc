/*
 *@BEGIN LICENSE
 *
 * cis by Psi4 Developer, a plugin to:
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
#include "psi4-dec.h"
#include <libdpd/dpd.h>
#include "psifiles.h"
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <boost/multi_array.hpp>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

INIT_PLUGIN

namespace psi{ namespace cis{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "CIS" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}


extern "C" PsiReturnType
cis(Options &options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");
   
    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int nirrep  = wfn->nirrep();
    int nmo     = wfn->nmo();
    //Assuming RHF reference
    int nocc = wfn->nalpha();
    int nsocc = 2 * nocc;
    int nso = 2 * nmo;
    int nvir = nso - nsocc;
    outfile->Printf("\nNocc:  %d  Nsocc: %d", nocc, nsocc);
    SharedVector motei(new Vector("MOTEI",nmo*nmo*nmo*nmo));
    SharedMatrix Hcis(new Matrix("Hcis", nsocc * nvir, nsocc*nvir));
    SharedMatrix F(new Matrix("F", nmo,nmo));
    SharedMatrix FSO(new Matrix("FSO", 2*nmo,2*nmo));
    SharedMatrix Ca(new Matrix("Ca", nmo,nmo));
    SharedMatrix H(new Matrix("H", nmo,nmo));
    
    typedef boost::multi_array<double, 4> array_type;
    
    array_type sotei(boost::extents[nso][nso][nso][nso]);
    array_type t2(boost::extents[nsocc][nsocc][nvir][nvir]);
    
    
    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");

    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                outfile->Printf( "(%2d %2d | %2d %2d) = %16.10f, "
                                 "symmetries = (%1d %1d | %1d %1d), "
                                 "relative indices = (%2d %2d | %2d %2d)\n",
                                 p, q, r, s, K.matrix[h][pq][rs], 
                                 psym, qsym, rsym, ssym, 
                                 prel, qrel, rrel, srel);
                int pqrs = INDEX4(p,q,r,s);
                motei->set(pqrs,K.matrix[h][pq][rs]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    double val1, val2;
    for(int p = 0; p < nso; p++){
       for(int q = 0; q < nso; q++){
          for(int r = 0; r < nso; r++){
             for(int s = 0; s < nso; s++){
                val1 = motei->get(INDEX4(p/2,r/2,q/2,s/2))*(p%2==r%2) * (q%2==s%2);
                val2 = motei->get(INDEX4(p/2,s/2,q/2,r/2))*(p%2==s%2)*(q%2==r%2);
                sotei[p][q][r][s] = val1 - val2;
             }
          }
       }
     }
    SharedMatrix Fmo(new Matrix("Fmo", nmo, nmo));        
    H= wfn->H();
    Fmo= wfn->Fa();
    Ca = wfn->Ca();
    H->transform(Ca);
    Fmo->transform(Ca);
    
    for(int i = 0; i < nso; i++){
       for(int j = 0; j < nso; j++){
          //FSO->set(i,j,H->get(i/2,j/2)*(i==j));
          FSO->set(i,j,Fmo->get(i/2,j/2)*(i==j));
          for(int m = 0; m < nsocc; m++){
          //   FSO->add(i,j,sotei[i][m][j][m]);
          }
       }
    }
    FSO->print();
    double val = 0.0;
    
    double mp2energy = 0.0;
    for(int i = 0; i < nsocc; i++){
       for(int a = 0; a < nvir; a++){
          for(int j = 0; j < nsocc; j++){
             for(int b = 0; b < nvir; b++){
                int ia = i*nvir + a;
                int jb = j*nvir + b;
                val = FSO->get(a+nsocc,b+nsocc)*(i==j) - FSO->get(i,j)*(a==b) + sotei[a+nsocc][j][i][b+nsocc];
                Hcis->set(ia,jb,val);
                outfile->Printf("\n\n ia = %d  jb = %d\t val = %20.12f", ia, jb,val);
                t2[i][j][a][b] = sotei[i][j][a+nsocc][b+nsocc]/(FSO->get(i,i) + FSO->get(j,j) - FSO->get(a + nsocc,a+nsocc) - FSO->get(b + nsocc, b+nsocc));
                outfile->Printf("\n\n\t t2[%d][%d][%d][%d] = %20.12f", i, j, a, b,t2[i][j][a][b]);
                mp2energy+=sotei[i][j][a+nsocc][b+nsocc]*(t2[i][j][a][b])*0.25;
              
              }
          }
       }
    }
    outfile->Printf("\n\n\t MP2 energy = %20.12f", mp2energy);
    Hcis->print();
    SharedMatrix Hevec(new Matrix("Hcis_evecs", nsocc*nvir, nsocc*nvir));
    SharedVector Heval(new Vector("HCis_evals", nsocc*nvir));
   
    Hcis->diagonalize(Hevec, Heval);
    Heval->print();
 
    return Success;
}

}} // End Namespaces
