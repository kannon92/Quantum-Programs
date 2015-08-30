/*
 *@BEGIN LICENSE
 *
 * scf by Psi4 Developer, a plugin to:
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
#include <libmints/vector.h>
#include <libpsio/psio.h>
#include "scf.h"
#include <cmath>

INIT_PLUGIN
namespace psi{ namespace scf {


extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "SCF"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
        options.add_bool("DO_FTHF", false);
        options.add_int("MAX_ITER", 50);
        options.add_int("MAX_TEMP",1000);
        options.add_int("TEMP_INCREMENT", 50);
    }

    return true;
}
Scfclass::Scfclass(boost::shared_ptr<Wavefunction> ref_wfn, Options& options) : Wavefunction(options, _default_psio_lib_)
{
    copy(ref_wfn);
}

extern "C"
PsiReturnType scf(Options &options)
{
     boost::shared_ptr<Scfclass> run_scf(new 
         Scfclass(Process::environment.wavefunction(), options));
     run_scf->ao_ints(options); 
     run_scf->scf_iteration();
     

    return Success;

}
Scfclass::~Scfclass()
{
}
void Scfclass::ao_ints(Options& options)
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    molecule->print();
    int nbf[] = { aoBasis->nbf() };
    double nucrep = molecule->nuclear_repulsion_energy();
    psi::outfile->Printf("\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);

    // The matrix factory can create matrices of the correct dimensions...
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, nbf, nbf);

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    sMat->print();
    sMat_= sMat->clone();
    tMat->print();
    vMat->print();

    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();
    hMat_ = hMat->clone();
    nbf_  = nbf[0];

    teiV_.reserve(nbf[0]*nbf[0]*nbf[0]*nbf[0]);
    if(doTei){

        // Now, the two-electron integrals
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        // The buffer will hold the integrals for each shell, as they're computed
        const double *buffer = eri->buffer();
        // The iterator conveniently lets us iterate over functions within shells
        AOShellCombinationsIterator shellIter = integral->shells_iterator();
        int count=0;
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // Compute quartet
            if (eri->compute_shell(shellIter)) {
                // From the quartet get all the integrals
                AOIntegralsIterator intIter = shellIter.integrals_iterator();
                for (intIter.first(); intIter.is_done() == false; intIter.next()) {
                    int p = intIter.i();
                    int q = intIter.j();
                    int r = intIter.k();
                    int s = intIter.l();
                    teiV_[INDEX4(p,q,r,s)]=buffer[intIter.index()];
                    

                }
            }
        }
        psi::outfile->Printf("\n\tThere are %d unique integrals\n\n", count);
    }
    

}
void Scfclass::scf_iteration()
{
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    int nocc = wfn->nalpha();
    nocc_ = nocc;
    outfile->Printf("\n nocc = %d", nocc);
    boost::shared_ptr<Matrix> Evec(new Matrix("Evectors", nbf_, nbf_));
    boost::shared_ptr<Vector> Eval(new Vector("Evals", nbf_));
   
    boost::shared_ptr<Matrix> Shalf(new Matrix("S^-(1/2)", nbf_, nbf_));
    boost::shared_ptr<Matrix> Svec(new Matrix("Evectors", nbf_, nbf_));
    boost::shared_ptr<Vector> Sval(new Vector("Evals", nbf_));
    boost::shared_ptr<Matrix> SI(new Matrix("SI", nbf_, nbf_));
    sMat_->diagonalize(Svec,Sval);
    for(int mu = 0; mu < nbf_; mu++){
       for(int nu = 0; nu < nbf_; nu++){
          int factor = mu==nu ? 1 : 0;
          SI->set(mu,nu,1.0/(sqrt(Sval->get(nu)))*factor);
       }
    }
    SI->print();
    for(int mu = 0; mu < nbf_; mu++){
       for(int nu = 0; nu < nbf_; nu++){
          for(int l = 0; l < nbf_; l++){
             for(int s = 0; s < nbf_; s++){
          Shalf->add(mu,nu,Svec->get(mu,l)*SI->get(l,s)*Svec->get(nu,s));
             }
          }
       } 
    } 
    Shalf->print();
    
    boost::shared_ptr<Matrix> Fold(new Matrix("Fock", nbf_, nbf_));
    Fold->copy(hMat_);
    Fold->transform(Shalf);
    Fold->print();
    
    SharedMatrix Fnotrans(new Matrix("Fnotrans", nbf_, nbf_)); 
    Fold->diagonalize(Evec,Eval);
 
    boost::shared_ptr<Matrix> C(new Matrix("C", nbf_, nbf_));
    boost::shared_ptr<Vector> eps(new Vector("epsilon", nbf_));

    C->gemm('F','F', 1.0, Shalf, Evec, 0.0);
    C->print();

    boost::shared_ptr<Matrix> D(new Matrix("D", nbf_,nbf_));
    boost::shared_ptr<Matrix> Dtest(new Matrix("Dtest", nbf_, nbf_));
   
    int iter = 0;
    bool t_run = false;
    if(options_.get_bool("DO_FTHF")==true){ 
      D = frac_occupation(C,Eval,iter, t_run);
    }
    else{
      D->gemm('N','T',nbf_,nbf_,nocc, 1.0,C, nbf_, C, nbf_, 0, nbf_, 0,0,0);
    }


    boost::shared_ptr<Matrix> FpH(new Matrix("FpH", nbf_,nbf_)); 
    FpH->copy(hMat_);
    FpH->add(hMat_);

    double energy_old = D->vector_dot(FpH);
    boost::shared_ptr<Molecule> mol = Process::environment.molecule();
    double nuc_rep = mol->nuclear_repulsion_energy();
    outfile->Printf("\nNuc Replusion = %20.12f", nuc_rep);
    outfile->Printf("\n 1st iteration = %20.12f", energy_old);
    int failiter = 0;
    double delta = energy_old;
    int maxiter = options_.get_int("MAX_ITER");
    //while((fabs(delta) > 1e-6) && failiter < maxiter)
 
    //This is used as a flag to state that T = 0 and ends iterations
    bool t_done = false;
    //Have not figured out why this loop doesn't work
    while(failiter < maxiter && std::fabs(delta) > 1e-8)
    {
        failiter++;
        
        energy_old = energy_;
        energy_ = 0.0;
       
        for(int mu = 0; mu < nbf_; mu++){
           for(int nu = 0; nu < nbf_; nu++){
              Fold->set(mu,nu,hMat_->get(mu,nu));
              for(int lam = 0; lam < nbf_; lam++){
                 for(int sig = 0; sig < nbf_; sig++){
                    Fold->add(mu,nu,D->get(lam,sig) * (2.0*teiV_[INDEX4(mu,nu,lam,sig)] - teiV_[INDEX4(mu,lam,nu,sig)]));
                 }
              }
           }
        }
        Fnotrans->copy(Fold);
        Fold->transform(Shalf); 
        Evec->zero();
        Eval->zero();
        Fold->diagonalize(Evec,Eval);
        C->zero();
        C->gemm(false,false, 1.0, Shalf, Evec,0.0);

      
        if(options_.get_bool("DO_FTHF")==true){
          D = frac_occupation(C,Eval, iter, t_done);
             if(t_done)
             {
                 outfile->Printf("\nI am here on %d\n", iter);
                 failiter = maxiter + 1;
             }
        }
        else
        {
             D->gemm('n','t',nbf_,nbf_,nocc, 1.0,C, nbf_, C, nbf_, 0, nbf_, 0,0,0);
        }
        //If the temperature is zero, do not calculation anything
        if(!t_done)
        {
           FpH->copy(hMat_);
           FpH->add(Fnotrans);
           energy_ = D->vector_dot(FpH);
        }
        else{
           outfile->Printf("\n Temperature is at 0.  I am taking a vacation \n");
           energy_ = D->vector_dot(FpH);
        }
   
        iter++;
        delta = energy_ - energy_old;  
        outfile->Printf("\n iter:%d  failiter:%d  energy:%20.12f  total:%20.12f  error:%20.12f", iter, failiter,energy_, energy_ + nuc_rep, energy_ - energy_old);
    }
    
   

}
double Scfclass::compute_energy()
{
   return 0.0;
}
boost::shared_ptr<Matrix> Scfclass::frac_occupation(SharedMatrix C, SharedVector e,int &iter, bool &tempdone)

{
    //P_{uv} = \sum n_i C C
    double T = 0.0;
    boost::shared_ptr<Matrix> D(new Matrix("DFT", nbf_, nbf_));

    //The maximum temperature.  High temperature should allow mixing between orbitals
    T = options_.get_int("MAX_TEMP");
    //The increment for decreasing temperature.  Borrowed from Scurseria fraction occupation convergence paper.  
    int increment = options_.get_int("TEMP_INCREMENT"); 
    T= T - iter*increment;

    //Conversion from K to atomic temperature.  
    T /=3.157746E5;
    std::vector<double> ni;

    //(HOMO + LUMO)/2 - nocc_ = number of orbitals occupied.  But, indexing starts from zero with SharedVector

    double ef = (e->get(nocc_ - 1) + e->get(nocc_))/2.0;

    outfile->Printf("\n Iter: %d  T: %12.3f  val:%20.6f",iter,T * 3.157746E5, ef);
    for(int i = 0; i < e->dim(); i++){
       //Fermi Dirac distribution - 1.0 / (1.0 + exp(\beta (e_i - ef))) 
       double val = 1.0/(1.0 + exp(1.0/(0.99994*T)*(e->get(i) - ef)));

       outfile->Printf("\n Orbital %d has an occupancy of %20.12f", i, val);
       ni.push_back(val);
    }
    for(int mu = 0; mu < nbf_; mu++){
       for(int nu = 0; nu < nbf_; nu++){
          for(int i = 0; i < nbf_; i++){
             D->add(mu,nu,ni[i]*C->get(mu,i)*C->get(nu,i)); 
          }
       }
    }
    
    if(T<= 0.00){
      //Need this to end
      tempdone= true;
    }
   
    return D;

}

}} // End Namespaces
