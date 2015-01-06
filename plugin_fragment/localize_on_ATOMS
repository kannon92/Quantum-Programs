#include <psi4-dec.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>

using namespace boost;

namespace psi{ namespace plugin_fragment {

void localize_on_atoms(Options& options, boost::shared_ptr<Wavefunction> wfn)
{
    int print = options.get_int("PRINT");

    // Compute the overlap matrix
    boost::shared_ptr<IntegralFactory> integral_ = wfn->integral();

    boost::shared_ptr<BasisSet> basisset_ = wfn->basisset();
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S_ao(new Matrix("S_ao",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S_ao);

    boost::shared_ptr<OneBodyAOInt> potential(integral_->ao_potential());
    SharedMatrix V_ao(new Matrix("V_ao",basisset_->nbf(),basisset_->nbf()));
    potential->compute(V_ao);

//    S_ao->print();
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    SharedMatrix AO2SO_ = pet->aotoso();
//    AO2SO_->print();

//    // Compute the W matrix for each fragment
//    std::vector<int> flist({0});
//    std::vector<int> glist;
//    boost::shared_ptr<Molecule> frag = mol->extract_subsets(flist,glist);
//    int frag_natom = frag->natom();
//    outfile->Printf("\n  Fragment contains %d atoms",frag->natom());

//    // Form a copy of S_ao and zero the rows and columns that are not on this fragment
////    max_a = min_a +
//    int nbfA = 0;

    std::vector<SharedMatrix> Ca_atoms;
    std::vector<SharedVector> lo_atoms;

    double threshold = 0.01;
    int natom = mol->natom();
    int nbf = basisset_->nbf();
    std::vector<int> atom_bf(natom,0);
    std::vector<int> atom_bf_begin(natom,0);
    std::vector<int> atom_bf_end(natom,0);
    int current_atom = -1;
    for (int mu = 0; mu < nbf; mu++) {
        int A = basisset_->function_to_center(mu);
        if (A > current_atom){
            atom_bf_begin[A] = mu;
            if (A >= 0){
                atom_bf_end[A-1] = mu;
            }
            current_atom = A;
        }
        atom_bf[A] += 1;
        outfile->Printf("\n  Function %d is on atom %d",mu,A);
    }
    atom_bf_end[natom-1] = nbf;

    for (int A = 0; A < natom; ++A){
        outfile->Printf("\n  Atom %2d : functions %2d  [%2d->%2d]",A,atom_bf[A],atom_bf_begin[A],atom_bf_end[A]-1);
    }

    // For now assume a closed-shell determinant
    int nocc = wfn->nalpha();
    int nvir = nbf - nocc;

    // Build a huge overlap matrix
    int nhuge = nocc * natom;
    SharedMatrix S_huge(new Matrix("S_huge",nhuge,nhuge));
    SharedMatrix U_huge(new Matrix("U_huge",nhuge,nhuge));
    SharedVector l_huge(new Vector("l_huge",nhuge));
    S_huge->zero();

    for (int A = 0; A < natom; ++A){

        outfile->Printf("\n\n  ==> Atom %d <==\n\n",A);
        // Form an atomic overlap matrix
        int nbfA = atom_bf[A];
        SharedMatrix S_f(new Matrix("S_f",nbfA,nbfA));
        for (int mu = 0; mu < nbfA; mu++) {
            for (int nu = 0; nu < nbfA; nu++) {
                S_f->set(mu,nu,S_ao->get(atom_bf_begin[A] + mu,atom_bf_begin[A] + nu));
            }
        }
//        S_f->print();

        // Invert S_f
        SharedMatrix S_f_inv(new Matrix("S_f_inv",nbfA,nbfA));
        S_f_inv->copy(S_f);
        S_f_inv->general_invert();
//        S_f_inv->print();

        // Place S_f_inv in a matrix of size nbf x
        SharedMatrix S_f_inv_large(new Matrix("S_f_inv_large",nbf,nbf));
        S_f_inv_large->zero();
        for (int mu = 0; mu < nbfA; mu++) {
            for (int nu = 0; nu < nbfA; nu++) {
                S_f_inv_large->set(atom_bf_begin[A] + mu,atom_bf_begin[A] + nu,S_f_inv->get(mu,nu));
            }
        }
//        S_f_inv_large->print();

        // Compute SA = C^T S S_f_inv S C
        SharedMatrix Ca = wfn->Ca();
        SharedMatrix SA(new Matrix("SA",nbf,nbf));
        SA->copy(S_f_inv_large);
        SA->transform(S_ao);
        SA->transform(Ca);
//        SA->print();

        // Extract the occupied part of SA and put it in S_huge
        SharedMatrix SAo(new Matrix("SAo",nocc,nocc));
        for (int i = 0; i < nocc; ++i){
            for (int j = 0; j < nocc; ++j){
                SAo->set(i,j,SA->get(i,j));
            }
        }
        SharedMatrix Uo(new Matrix("Uo",nocc,nocc));
        SharedVector lo(new Vector("lo",nocc));
        SAo->diagonalize(Uo,lo,descending);
//        lo->print();

        SharedMatrix SAv(new Matrix("SAv",nvir,nvir));
        for (int i = 0; i < nvir; ++i){
            for (int j = 0; j < nvir; ++j){
                SAv->set(i,j,SA->get(i + nocc,j + nocc));
            }
        }
        SharedMatrix Uv(new Matrix("Uv",nvir,nvir));
        SharedVector lv(new Vector("lv",nvir));
        SAv->diagonalize(Uv,lv,descending);
//        lv->print();

        SharedMatrix U(new Matrix("U",nbf,nbf));
        for (int i = 0; i < nocc; ++i){
            for (int j = 0; j < nocc; ++j){
                U->set(i,j,Uo->get(i,j));
            }
        }
        for (int a = 0; a < nvir; ++a){
            for (int b = 0; b < nvir; ++b){
                U->set(a + nocc,b + nocc,Uv->get(a,b));
            }
        }
//        U->print();
        SharedMatrix Ca_new(Ca->clone());
        Ca_new->gemm(false, false, 1.0,Ca,U, 0.0);
        Ca_atoms.push_back(Ca_new);
        lo_atoms.push_back(lo);
    }

//    int nloc_occ = 0;
//    std::vector<std::pair<int,int>> loc_occ;
//    for (int A = 0; A < natom; ++A){
//        for (int i = 0; i < nocc; ++i){
//            double lo = lo_atoms[A]->get(i);
//            if (std::fabs(lo) > threshold){
//                outfile->Printf("\n  Adding atomic localize occupied orbital %2d from atom %2d (%f)",i,A,lo);
//                nloc_occ += 1;
//                loc_occ.push_back(std::make_pair(A,i));
//            }
//        }
//    }
//    outfile->Printf("\n  Added a total of %d orbitals",nloc_occ);


    // Sort the orbitals according to the value of lambda_o
    std::vector<std::tuple<double,int,int>> sorted_loc_occ;
    for (int A = 0; A < natom; ++A){
        for (int i = 0; i < nocc; ++i){
            double lo = lo_atoms[A]->get(i);
            sorted_loc_occ.push_back(std::make_tuple(lo,A,i));
        }
    }
    std::sort(sorted_loc_occ.begin(),sorted_loc_occ.end());
    std::reverse(sorted_loc_occ.begin(),sorted_loc_occ.end());

    // Select the orbitals to keep
    int nkeep = 0;
    for (int t = 0; t < natom * nocc; ++t){
        double lo;
        int A,i;
        // access the elements by unpacking them with std::tie ...
        std::tie(lo,A,i) = sorted_loc_occ[t];
        outfile->Printf("\n %3d  %2d-%2d : %f",t+1,A,i,lo);
        if (lo > threshold){
            nkeep += 1;
            outfile->Printf(" --> keep");
        }else{
        }
    }
    outfile->Printf("\n Keeping %d MOs",nkeep);

    // Place the orbitals to keep in Cloc
    SharedMatrix Cloc(new Matrix("Cloc",nbf,nkeep));
    for (int t = 0; t < nkeep; ++t){
        double lo;
        int A,i;
        // access the elements by unpacking them with std::tie ...
        std::tie(lo,A,i) = sorted_loc_occ[t];
        for (int mu = 0; mu < nbf; ++mu){
            Cloc->set(mu,t,Ca_atoms[A]->get(mu,i));
        }
    }
//    Cloc->print();


    SharedMatrix Ca = wfn->Ca();

    int procedure = 1;
    // Form the overlap matrix and fiddle with it
    if(procedure == 0){
        SharedMatrix XX(S_ao->clone());
        XX->transform(Cloc);

        // Modify the diagonal
        for (int t = 0; t < nkeep; ++t){
            double lo;
            int A,i;
            // access the elements by unpacking them with std::tie ...
            std::tie(lo,A,i) = sorted_loc_occ[t];
            double diag = XX->get(t,t);
//            XX->set(t,t,diag / lo);
        }
//        XX->print();

        // Diagonalize XX
        SharedMatrix U(new Matrix("OmegaU",nkeep,nkeep));
        SharedVector lambda(new Vector("Omegal",nkeep));
        XX->diagonalize(U,lambda,descending);
//        lambda->print();

        SharedMatrix X(new Matrix("X",nkeep,nocc));
        for (int i = 0; i < nkeep; ++i){
            for (int j = 0; j < nocc; ++j){
                X->set(i,j,U->get(i,j)/std::sqrt(lambda->get(j)));
            }
        }
        SharedMatrix Cloc_ortho(new Matrix("Cloc_ortho",nbf,nocc));
        Cloc_ortho->gemm(false, false, 1.0,Cloc,X, 0.0);

        SharedMatrix XX2(S_ao->clone());
        XX2->transform(Cloc_ortho);
//        XX2->print();
    }

    // Gram-Schmidt
    if(procedure == 1){
        outfile->Printf("\n  Performing Schmidt orthogonalization");
//        SharedMatrix XX(S_ao->clone());
        SharedMatrix XX(V_ao->clone());
        XX->transform(Cloc);
//        XX->print();

        SharedMatrix Cloc_ortho(new Matrix("Cloc_ortho",nbf,nocc));
        double gs_threshold = 1.0e-9;
        int nonzero = 0;
        for (int i = 0; i < nkeep; ++i){
            // Compute the norm of all orbitals
            XX->copy(V_ao);
            XX->transform(Cloc);

            // Find the orbital with the largest norm
            std::vector<std::pair<double,int>> sorted_gs;
            for (int j = i; j < nkeep; ++j){
                sorted_gs.push_back(std::make_pair(XX->get(j,j),j));
                outfile->Printf("\n  XX %d -> %f",j,XX->get(j,j));
            }
            std::sort(sorted_gs.begin(),sorted_gs.end());
//            std::reverse(sorted_gs.begin(),sorted_gs.end());

            int ii = sorted_gs[0].second;
            double norm = 0.0;
            for (int mu = 0; mu < nbf; ++mu){
                for (int nu = 0; nu < nbf; ++nu){
                    norm += Cloc->get(mu,ii) * S_ao->get(mu,nu) * Cloc->get(nu,ii);
                }
            }
            norm = std::sqrt(norm);

            outfile->Printf("\n  GS %d -> %d (norm = %f)",i,ii,norm);
            if (norm > gs_threshold){
                // Normalize
                for (int mu = 0; mu < nbf; ++mu){
                    Cloc_ortho->set(mu,nonzero,Cloc->get(mu,ii) / norm);
                }
                // Project out from the others
                for (int j = 0; j < nkeep; ++j){
                    double proj = 0.0;
                    for (int mu = 0; mu < nbf; ++mu){
                        for (int nu = 0; nu < nbf; ++nu){
                            proj += Cloc_ortho->get(mu,nonzero) * S_ao->get(mu,nu) * Cloc->get(nu,j);
                        }
                    }
                    for (int mu = 0; mu < nbf; ++mu){
                        Cloc->set(mu,j,Cloc->get(mu,j) - proj * Cloc_ortho->get(mu,nonzero));
                    }
                }
                nonzero += 1;
            }
            if(nonzero == nocc) break;
        }

        for (int mu = 0; mu < nbf; ++mu){
            for (int i = 0; i < nocc; ++i){
                Ca->set(mu,i,Cloc_ortho->get(mu,i));
            }
        }
        SharedMatrix XX2(S_ao->clone());
        XX2->transform(Cloc_ortho);
//        XX2->print();
    }



//    for (int mu = 0; mu < nbf; ++mu){
//        for (int i = 0; i < nocc; ++i){
//            Ca->set(mu,i,Cloc->get(mu,i));
//            //Ca->set(mu,i,Cloc_ortho->get(mu,i));
////            Ca->set(mu,i,Ca_atoms->ge);
//        }
//    }


//    Procedure for orthogonalization
//    SharedMatrix Cnu(new Matrix("Cnu",nkeep,nkeep));
//    SharedMatrix M(new Matrix("M",nkeep,nkeep));
//    SharedMatrix Omega(new Matrix("Omega",nkeep,nkeep));
//    SharedMatrix A(new Matrix("A",nkeep,nkeep));

//    Cnu->zero();
//    for (int i = 0; i < nkeep; ++i){
//        Cnu->set(i,i,1.0);
//    }

//    SharedMatrix XX(S_ao->clone());
//    XX->transform(Cloc);
//    XX->print();

//    std::vector<double> oldC(nkeep);
//    for (int cycle = 0; cycle < 10; ++cycle){
//        M->copy(XX);
//        M->transform(Cnu);
////        M->print();
//        Omega->copy(M);
//        Omega->power(0.5);
////        SharedMatrix OmegaU(new Matrix("OmegaU",nkeep,nkeep));
////        SharedVector Omegal(new Vector("Omegal",nkeep));
////        XX->diagonalize(OmegaU,Omegal);
////        Omegal->print();
////        Omega->print();
////        Cnu->print();
//        double rms = 0.0;
//        for (int i = 0; i < nkeep; ++i){
//            double c = Omega->get(i,i) / Cnu->get(i,i);
//            Cnu->set(i,i,c);
//            rms += std::pow(c-oldC[i],2.0);
//            oldC[i] = c;
//        }
//        outfile->Printf("\n %d -> %20.12f",cycle,std::sqrt(rms));
//    }
//    SharedMatrix Omega_inv(Omega->clone());
//    Omega_inv->pseudoinverse(1.0e-6);
//    Omega_inv->print();


//    A->gemm(false, false, 1.0,Omega,Omega_inv, 0.0);
//    A->print();
//    A->gemm(false, false, 1.0,Cnu,Omega_inv, 0.0);

//    SharedMatrix Cloc_ortho(new Matrix("Cloc_ortho",nbf,nkeep));
//    Cloc_ortho->gemm(false, false, 1.0,Cloc,A, 0.0);

//    SharedMatrix test(S_ao->clone());
//    test->transform(Cloc_ortho);
//    test->print();

}

}} // End namespaces
