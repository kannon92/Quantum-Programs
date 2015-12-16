/*
 *@BEGIN LICENSE
 *
 * integral_trans by Psi4 Developer, a plugin to:
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
#include <libfock/jk.h>
#include <libmints/mints.h>
#include <libqt/qt.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

INIT_PLUGIN

namespace psi{ namespace integral_trans{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "INTEGRAL_TRANS" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}


extern "C" PsiReturnType
integral_trans(Options &options)
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

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    SharedMatrix TEI_PSI4(new Matrix("TEI", wfn->nmo() * wfn->nmo(), wfn->nmo() * wfn->nmo()));
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
                psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f, "
                                 "symmetries = (%1d %1d | %1d %1d), "
                                 "relative indices = (%2d %2d | %2d %2d)\n",
                                 p, q, r, s, K.matrix[h][pq][rs], 
                                 psym, qsym, rsym, ssym, 
                                 prel, qrel, rrel, srel);
                TEI_PSI4->set(pq, rs, K.matrix[h][pq][rs]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    outfile->Printf("\n TEI_PSI4: %8.8f", TEI_PSI4->rms());

    /// Do the MO transform using JK builders

    SharedMatrix C = wfn->Ca();
    std::vector<SharedMatrix> D_vec;
    int nmo = wfn->nmo();
    SharedMatrix Identity(new Matrix("I", nmo, nmo));
    Identity->identity();

    std::vector<SharedMatrix> J_pq;
    SharedMatrix TEI_MO(new Matrix("TEI_MO", nmo * nmo, nmo * nmo));
    for(int i = 0; i < nmo; i++){
        SharedVector C_i = C->get_column(0, i);
        for(int j = 0; j < nmo; j++){
            SharedMatrix D(new Matrix("D", nmo, nmo));
            SharedVector C_j = C->get_column(0, j);
            //C_j->print();
            C_DGER(nmo, nmo, 1.0, &(C_i->pointer()[0]), 1, &(C_j->pointer()[0]), 1, D->pointer()[0], nmo);
            ///Form D = C_rC_s'
            ///Use this to compute J matrix
            /// (pq | rs) = C J C'

            boost::shared_ptr<JK> JK_core = JK::build_JK();
            JK_core->set_memory(Process::environment.get_memory() * 0.8);
            JK_core->initialize();

            std::vector<boost::shared_ptr<Matrix> > &Cl = JK_core->C_left();
            std::vector<boost::shared_ptr<Matrix> > &Cr = JK_core->C_right();
            Cl.push_back(D);
            Cr.push_back(Identity);
            JK_core->set_do_K(false);
            JK_core->compute();

            SharedMatrix J = JK_core->J()[0];
            J->transform(C);
            for(int p = 0; p < nmo; p++){
                for(int q = 0; q < nmo; q++){
                    TEI_MO->set(p * nmo + q, i * nmo + j, J->get(p, q));
                }
            }
        }
    }
    outfile->Printf("\n\n TEI_MO norm: %8.8f", TEI_MO->rms());


    return Success;
}

}} // End Namespaces
