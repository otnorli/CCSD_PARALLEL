/*
   This file is part of toffyrn::scotch.

   toffyrn::scotch is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   toffyrn::scotch is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with toffyrn::scotch.  If not, see <http://www.gnu.org/licenses/>.
 */
/* 
 * File:   HFbasis.cpp
 * Author: toffyrn
 * 
 * Created on 16. februar 2012, 17:21
 */

#include "HFsys.h"
#include "System.h"
#include "GEMM.h"
#include "CLgemm.h"

using namespace toffyrn::libUNK;
using namespace arma;
using namespace std;

HFsys::HFsys(System const * originalSys, arma::mat const &Coeff, GEMM * mult)
{
    if (mult == NULL) //Use standard GEMM object.
        this->mult = &mult_standard; 
    else
        this->mult = mult;

    this->basis = originalSys->get_basis()->clone();
    this->f_hh = *originalSys->get_f_hh();
    this->f_ph = *originalSys->get_f_ph();
    this->f_pp = *originalSys->get_f_pp();
    this->v_hhhh = *originalSys->get_v_hhhh();
    this->v_phhh = *originalSys->get_v_phhh();
    this->v_pphh = *originalSys->get_v_pphh();
    this->v_phph = *originalSys->get_v_phph();
    this->v_ppph = *originalSys->get_v_ppph();
    this->v_pppp = *originalSys->get_v_pppp();
    transformElements(Coeff, originalSys);
}

HFsys::~HFsys()
{

}

double HFsys::f_elem(std::size_t p, std::size_t q) const
{
    int nH = basis->get_nH();

    if (p < nH && q < nH)
        return f_hh(p, q);
    else if (p >= nH && q < nH)
        return f_ph(p - nH, q);
    else if (p < nH && q >= nH)
        return f_ph(q - nH, p);
    else if (p >= nH && q >= nH)
        return f_pp(p - nH, q - nH);
    else
        throw std::string("Unreachable element in f_elem.");
}

double HFsys::v_elem(
        std::size_t p, std::size_t q,
        std::size_t r, std::size_t s) const
{
    //Needed mappings & constants
    int nH = basis->get_nH();
    int nP = basis->get_nP();
    umat const * mapHH = basis->get_map_lm_lmdMU();
    umat const * mapPH = basis->get_map_dl_lmdNU();
    umat const * mapPP = basis->get_map_de_lmdXI();

    //Element value and phasefactor
    double elem = 0;
    double phase = 1.0;

    //Is there 0, 1 or 2 particle states in bra or ket?
    int pStatesLeft = 0;
    if (p >= nH)
        pStatesLeft++;
    if (q >= nH)
        pStatesLeft++;
    int pStatesRight = 0;
    if (r >= nH)
        pStatesRight++;
    if (s >= nH)
        pStatesRight++;

    //Gather particle states to the left.
    if (pStatesRight > pStatesLeft)
    {
        //Swap p,q <-> r,s
        int tmp = p;
        p = r;
        r = tmp;
        tmp = q;
        q = s;
        s = tmp;
        //Left <-> Right
        tmp = pStatesRight;
        pStatesRight = pStatesLeft;
        pStatesLeft = tmp;
    }
    if (q > p)
    {
        int tmp = p;
        p = q;
        q = tmp;
        phase *= -1;
    }
    if (s > r)
    {
        int tmp = r;
        r = s;
        s = tmp;
        phase *= -1;
    }



    if (pStatesLeft == 0 && pStatesRight == 0)
    { //hhhh
        int pq = p + q * nH;
        int lmd_pq = (*mapHH)(0, pq);
        int mu_pq = (*mapHH)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_hhhh.at(lmd_pq)(mu_pq, mu_rs);
    } else if (pStatesLeft == 1 && pStatesRight == 0)
    { //phhh
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*mapPH)(0, pq);
        int nu_pq = (*mapPH)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_phhh.at(lmd_pq)(nu_pq, mu_rs);
    } else if (pStatesLeft == 1 && pStatesRight == 1)
    { //phph
        int pq = (p - nH) + q * nP;
        int lmd_pq = (*mapPH)(0, pq);
        int nu_pq = (*mapPH)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*mapPH)(0, rs);
        int nu_rs = (*mapPH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_phph.at(lmd_pq)(nu_pq, nu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 0)
    { //pphh
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = r + s * nH;
        int lmd_rs = (*mapHH)(0, rs);
        int mu_rs = (*mapHH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_pphh.at(lmd_pq)(xi_pq, mu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 1)
    { //ppph
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = (r - nH) + s * nP;
        int lmd_rs = (*mapPH)(0, rs);
        int nu_rs = (*mapPH)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_ppph.at(lmd_pq)(xi_pq, nu_rs);
    } else if (pStatesLeft == 2 && pStatesRight == 2)
    { //pppp
        int pq = (p - nH) + (q - nH) * nP;
        int lmd_pq = (*mapPP)(0, pq);
        int xi_pq = (*mapPP)(1, pq);
        int rs = (r - nH) + (s - nH) * nP;
        int lmd_rs = (*mapPP)(0, rs);
        int xi_rs = (*mapPP)(1, rs);
        if (lmd_pq == lmd_rs)
            elem = v_pppp.at(lmd_pq)(xi_pq, xi_rs);
    }

    return phase * elem;
}

void HFsys::transformElements(mat const &C, System const * origSys)
{
    wall_clock timer;
    timer.tic();

    int nH = basis->get_nH();
    int nP = basis->get_nP();
    int nTot = nH + nP;
    int dimLMD2 = basis->dim_lmd_2p();
    vector<uvec> const * mapHH = basis->get_map_lmdMU_lm();
    vector<uvec> const * mapPH = basis->get_map_lmdNU_dl();
    vector<uvec> const * mapPP = basis->get_map_lmdXI_de();

    //Get the singleparticle interactions.
    mat h0(nTot, nTot);
    span hSpan(0, nH - 1);
    span pSpan(nH, nTot - 1);
    h0(hSpan, hSpan) = f_hh;
    h0(pSpan, hSpan) = f_ph;
    h0(hSpan, pSpan) = trans(h0(pSpan, hSpan));
    h0(pSpan, pSpan) = f_pp;
    for (int p = 0; p < nTot; p++)
        for (int q = 0; q < nTot; q++)
        {
            double u_pq = 0;
            for (int i = 0; i < nH; i++)
                u_pq += v_elem(p, i, q, i);
            h0(p, q) -= u_pq;
        }
    //Transform h0
    h0 = C * h0 * C.t();

    //Old matrices 
    vector<mat> const * old_hhhh = origSys->get_v_hhhh();
    vector<mat> const * old_phhh = origSys->get_v_phhh();
    vector<mat> const * old_pphh = origSys->get_v_pphh();
    vector<mat> const * old_phph = origSys->get_v_phph();
    vector<mat> const * old_ppph = origSys->get_v_ppph();
    vector<mat> const * old_pppp = origSys->get_v_pppp();

    //Transform one channel at a time
    for (int lmd2 = 0; lmd2 < dimLMD2; lmd2++)
    {
        int dimMU = mapHH->at(lmd2).size();
        int dimNU = mapPH->at(lmd2).size();
        int dimXI = mapPP->at(lmd2).size();

        //Fill
        mat C_hhhh = zeros<mat > (dimMU, dimMU);
        for (int mu1 = 0; mu1 < dimMU; mu1++)
            for (int mu2 = 0; mu2 < dimMU; mu2++)
            {
                int ij = mapHH->at(lmd2)(mu1);
                int i = ij % nH;
                int j = ij / nH;
                int kl = mapHH->at(lmd2)(mu2);
                int k = kl % nH;
                int l = kl / nH;
                C_hhhh(mu1, mu2) = C(i, k) * C(j, l);
            }
        mat C_hhph = zeros<mat > (dimMU, dimNU);
        for (int mu = 0; mu < dimMU; mu++)
            for (int nu = 0; nu < dimNU; nu++)
            {
                int ij = mapHH->at(lmd2)(mu);
                int i = ij % nH;
                int j = ij / nH;
                int cl = mapPH->at(lmd2)(nu);
                int c = cl % nP;
                int l = cl / nP;
                C_hhph(mu, nu) = C(i, c + nH) * C(j, l);
            }
        mat C_hhhp = zeros<mat > (dimMU, dimNU);
        for (int mu = 0; mu < dimMU; mu++)
            for (int nu = 0; nu < dimNU; nu++)
            {
                int ij = mapHH->at(lmd2)(mu);
                int i = ij % nH;
                int j = ij / nH;
                int dk = mapPH->at(lmd2)(nu);
                int d = dk % nP;
                int k = dk / nP;
                C_hhhp(mu, nu) = C(i, k) * C(j, d + nH);
            }
        mat C_phhh = zeros<mat > (dimNU, dimMU);
        for (int nu = 0; nu < dimNU; nu++)
            for (int mu = 0; mu < dimMU; mu++)
            {
                int aj = mapPH->at(lmd2)(nu);
                int a = aj % nP;
                int j = aj / nP;
                int kl = mapHH->at(lmd2)(mu);
                int k = kl % nH;
                int l = kl / nH;
                C_phhh(nu, mu) = C(a + nH, k) * C(j, l);
            }
        mat C_phhp = zeros<mat > (dimNU, dimNU);
        for (int nu1 = 0; nu1 < dimNU; nu1++)
            for (int nu2 = 0; nu2 < dimNU; nu2++)
            {
                int aj = mapPH->at(lmd2)(nu1);
                int a = aj % nP;
                int j = aj / nP;
                int dk = mapPH->at(lmd2)(nu2);
                int d = dk % nP;
                int k = dk / nP;
                C_phhp(nu1, nu2) = C(a + nH, k) * C(j, d + nH);
            }
        mat C_phpp = zeros<mat > (dimNU, dimXI);
        for (int nu = 0; nu < dimNU; nu++)
            for (int xi = 0; xi < dimXI; xi++)
            {
                int aj = mapPH->at(lmd2)(nu);
                int a = aj % nP;
                int j = aj / nP;
                int cd = mapPP->at(lmd2)(xi);
                int c = cd % nP;
                int d = cd / nP;
                C_phpp(nu, xi) = C(a + nH, c + nH) * C(j, d + nH);
            }
        mat C_pphh = zeros<mat > (dimXI, dimMU);
        for (int xi = 0; xi < dimXI; xi++)
            for (int mu = 0; mu < dimMU; mu++)
            {
                int ab = mapPP->at(lmd2)(xi);
                int a = ab % nP;
                int b = ab / nP;
                int kl = mapHH->at(lmd2)(mu);
                int k = kl % nH;
                int l = kl / nH;
                C_pphh(xi, mu) = C(a + nH, k) * C(b + nH, l);
            }
        mat C_pphp = zeros<mat > (dimXI, dimNU);
        for (int xi = 0; xi < dimXI; xi++)
            for (int nu = 0; nu < dimNU; nu++)
            {
                int ab = mapPP->at(lmd2)(xi);
                int a = ab % nP;
                int b = ab / nP;
                int dk = mapPH->at(lmd2)(nu);
                int d = dk % nP;
                int k = dk / nP;
                C_pphp(xi, nu) = C(a + nH, k) * C(b + nH, d + nH);
            }
        mat C_hhpp = zeros<mat > (dimMU, dimXI);
        for (int mu = 0; mu < dimMU; mu++)
            for (int xi = 0; xi < dimXI; xi++)
            {
                int ij = mapHH->at(lmd2)(mu);
                int i = ij % nH;
                int j = ij / nH;
                int cd = mapPP->at(lmd2)(xi);
                int c = cd % nP;
                int d = cd / nP;
                C_hhpp(mu, xi) = C(i, c + nH) * C(j, d + nH);
            }
        mat C_phph = zeros<mat > (dimNU, dimNU);
        for (int nu1 = 0; nu1 < dimNU; nu1++)
            for (int nu2 = 0; nu2 < dimNU; nu2++)
            {
                int aj = mapPH->at(lmd2)(nu1);
                int a = aj % nP;
                int j = aj / nP;
                int cl = mapPH->at(lmd2)(nu2);
                int c = cl % nP;
                int l = cl / nP;
                C_phph(nu1, nu2) = C(a + nH, c + nH) * C(j, l);
            }
        mat C_ppph = zeros<mat > (dimXI, dimNU);
        for (int xi = 0; xi < dimXI; xi++)
            for (int nu = 0; nu < dimNU; nu++)
            {
                int ab = mapPP->at(lmd2)(xi);
                int a = ab % nP;
                int b = ab / nP;
                int cl = mapPH->at(lmd2)(nu);
                int c = cl % nP;
                int l = cl / nP;
                C_ppph(xi, nu) = C(a + nH, c + nH) * C(b + nH, l);
            }
        mat C_pppp = zeros<mat > (dimXI, dimXI);
        for (int xi1 = 0; xi1 < dimXI; xi1++)
            for (int xi2 = 0; xi2 < dimXI; xi2++)
            {
                int ab = mapPP->at(lmd2)(xi1);
                int a = ab % nP;
                int b = ab / nP;
                int cd = mapPP->at(lmd2)(xi2);
                int c = cd % nP;
                int d = cd / nP;
                C_pppp(xi1, xi2) = C(a + nH, c + nH) * C(b + nH, d + nH);
            }

        //Intermediate multiplications
        mat int_Cpppp_Vpppp = mult->dgemm(C_pppp, old_pppp->at(lmd2));
        mat int_Cpphp_Vpppht = mult->dgemm(C_pphp, old_ppph->at(lmd2).t());
        mat int_Cppph_Vpppht = mult->dgemm(C_ppph, old_ppph->at(lmd2).t());
        mat int_Cpppp_Vppph = mult->dgemm(C_pppp, old_ppph->at(lmd2));
        mat int_Cpphh_Vhhhh = mult->dgemm(C_pphh, old_hhhh->at(lmd2));
        mat int_Cppph_Vphhh = mult->dgemm(C_ppph, old_phhh->at(lmd2));
        mat int_Cpphp_Vphhh = mult->dgemm(C_pphp, old_phhh->at(lmd2));
        mat int_Cpphh_Vphhht = mult->dgemm(C_pphh, old_phhh->at(lmd2).t());
        mat int_Cpppp_Vpphh = mult->dgemm(C_pppp, old_pphh->at(lmd2));
        mat int_Cpphh_Vpphht = mult->dgemm(C_pphh, old_pphh->at(lmd2).t());
        mat int_Cppph_Vphph = mult->dgemm(C_ppph, old_phph->at(lmd2));
        mat int_Cpphp_Vphph = mult->dgemm(C_pphp, old_phph->at(lmd2));
        mat int_Cphhh_Vhhhh = mult->dgemm(C_phhh, old_hhhh->at(lmd2));
        mat int_Cphph_Vphhh = mult->dgemm(C_phph, old_phhh->at(lmd2));
        mat int_Cphhp_Vphhh = mult->dgemm(C_phhp, old_phhh->at(lmd2));
        mat int_Cphhh_Vphhht = mult->dgemm(C_phhh, old_phhh->at(lmd2).t());
        mat int_Cphpp_Vpphh = mult->dgemm(C_phpp, old_pphh->at(lmd2));
        mat int_Cphhh_Vpphht = mult->dgemm(C_phhh, old_pphh->at(lmd2).t());
        mat int_Cphph_Vphph = mult->dgemm(C_phph, old_phph->at(lmd2));
        mat int_Cphhp_Vphph = mult->dgemm(C_phhp, old_phph->at(lmd2));
        mat int_Cphpp_Vppph = mult->dgemm(C_phpp, old_ppph->at(lmd2));
        mat int_Cphph_Vpppht = mult->dgemm(C_phph, old_ppph->at(lmd2).t());
        mat int_Cphhp_Vpppht = mult->dgemm(C_phhp, old_ppph->at(lmd2).t());
        mat int_Cphpp_Vpppp = mult->dgemm(C_phpp, old_pppp->at(lmd2));
        mat int_Chhhh_Vhhhh = mult->dgemm(C_hhhh, old_hhhh->at(lmd2));
        mat int_Chhph_Vphhh = mult->dgemm(C_hhph, old_phhh->at(lmd2));
        mat int_Chhhp_Vphhh = mult->dgemm(C_hhhp, old_phhh->at(lmd2));
        mat int_Chhhh_Vphhht = mult->dgemm(C_hhhh, old_phhh->at(lmd2).t());
        mat int_Chhpp_Vpphh = mult->dgemm(C_hhpp, old_pphh->at(lmd2));
        mat int_Chhhh_Vpphht = mult->dgemm(C_hhhh, old_pphh->at(lmd2).t());
        mat int_Chhph_Vphph = mult->dgemm(C_hhph, old_phph->at(lmd2));
        mat int_Chhhp_Vphph = mult->dgemm(C_hhhp, old_phph->at(lmd2));
        mat int_Chhpp_Vppph = mult->dgemm(C_hhpp, old_ppph->at(lmd2));
        mat int_Chhph_Vpppht = mult->dgemm(C_hhph, old_ppph->at(lmd2).t());
        mat int_Chhhp_Vpppht = mult->dgemm(C_hhhp, old_ppph->at(lmd2).t());
        mat int_Chhpp_Vpppp = mult->dgemm(C_hhpp, old_pppp->at(lmd2));

        //hhhh
        {
            //0p4h
            v_hhhh.at(lmd2) = mult->dgemm(int_Chhhh_Vhhhh, C_hhhh.t());
            //1p3h
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhph_Vphhh, C_hhhh.t());
            v_hhhh.at(lmd2) -= mult->dgemm(int_Chhhp_Vphhh, C_hhhh.t());
            //            new_hhhh.at(lmd2) += int_Chhhh_Vphhht * C_hhph.t();
            //            new_hhhh.at(lmd2) -= int_Chhhh_Vphhht * C_hhhp.t();
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhhh_Vphhht, (C_hhph.t() - C_hhhp.t()));
            //2p2h
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhpp_Vpphh, C_hhhh.t());
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhhh_Vpphht, C_hhpp.t());
            //            new_hhhh.at(lmd2) += int_Chhph_Vphph * C_hhph.t();
            //            new_hhhh.at(lmd2) -= int_Chhph_Vphph * C_hhhp.t();
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhph_Vphph, (C_hhph.t() - C_hhhp.t()));
            //            new_hhhh.at(lmd2) -= int_Chhhp_Vphph * C_hhph.t();
            //            new_hhhh.at(lmd2) += int_Chhhp_Vphph * C_hhhp.t();
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhhp_Vphph, (C_hhhp.t() - C_hhph.t()));
            //3p1h
            //            new_hhhh.at(lmd2) += int_Chhpp_Vppph * C_hhph.t();
            //            new_hhhh.at(lmd2) -= int_Chhpp_Vppph * C_hhhp.t();
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhpp_Vppph, (C_hhph.t() - C_hhhp.t()));
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhph_Vpppht, C_hhpp.t());
            v_hhhh.at(lmd2) -= mult->dgemm(int_Chhhp_Vpppht, C_hhpp.t());
            //4p0h
            v_hhhh.at(lmd2) += mult->dgemm(int_Chhpp_Vpppp, C_hhpp.t());
        }
        //phhh
        {
            //0p4h
            v_phhh.at(lmd2) = mult->dgemm(int_Cphhh_Vhhhh, C_hhhh.t());
            //1p3h
            v_phhh.at(lmd2) += mult->dgemm(int_Cphph_Vphhh, C_hhhh.t());
            v_phhh.at(lmd2) -= mult->dgemm(int_Cphhp_Vphhh, C_hhhh.t());
            //            new_phhh.at(lmd2) += int_Cphhh_Vphhht * C_hhph.t();
            //            new_phhh.at(lmd2) -= int_Cphhh_Vphhht * C_hhhp.t();
            v_phhh.at(lmd2) += mult->dgemm(int_Cphhh_Vphhht, (C_hhph.t() - C_hhhp.t()));
            //2p2h
            v_phhh.at(lmd2) += mult->dgemm(int_Cphpp_Vpphh, C_hhhh.t());
            v_phhh.at(lmd2) += mult->dgemm(int_Cphhh_Vpphht, C_hhpp.t());
            //            new_phhh.at(lmd2) += int_Cphph_Vphph * C_hhph.t();
            //            new_phhh.at(lmd2) -= int_Cphph_Vphph * C_hhhp.t();
            v_phhh.at(lmd2) += mult->dgemm(int_Cphph_Vphph, (C_hhph.t() - C_hhhp.t()));
            //            new_phhh.at(lmd2) -= int_Cphhp_Vphph * C_hhph.t();
            //            new_phhh.at(lmd2) += int_Cphhp_Vphph * C_hhhp.t();
            v_phhh.at(lmd2) += mult->dgemm(int_Cphhp_Vphph, (C_hhhp.t() - C_hhph.t()));
            //3p1h
            //            new_phhh.at(lmd2) += int_Cphpp_Vppph * C_hhph.t();
            //            new_phhh.at(lmd2) -= int_Cphpp_Vppph * C_hhhp.t();
            v_phhh.at(lmd2) += mult->dgemm(int_Cphpp_Vppph, (C_hhph.t() - C_hhhp.t()));
            v_phhh.at(lmd2) += mult->dgemm(int_Cphph_Vpppht, C_hhpp.t());
            v_phhh.at(lmd2) -= mult->dgemm(int_Cphhp_Vpppht, C_hhpp.t());
            //4p0h
            v_phhh.at(lmd2) += mult->dgemm(int_Cphpp_Vpppp, C_hhpp.t());
        }
        //pphh
        {
            //0p4h
            v_pphh.at(lmd2) = mult->dgemm(int_Cpphh_Vhhhh, C_hhhh.t());
            //1p3h
            v_pphh.at(lmd2) += mult->dgemm(int_Cppph_Vphhh, C_hhhh.t());
            v_pphh.at(lmd2) -= mult->dgemm(int_Cpphp_Vphhh, C_hhhh.t());
            //            new_pphh.at(lmd2) += int_Cpphh_Vphhht * C_hhph.t();
            //            new_pphh.at(lmd2) -= int_Cpphh_Vphhht * C_hhhp.t();
            v_pphh.at(lmd2) += mult->dgemm(int_Cpphh_Vphhht, (C_hhph.t() - C_hhhp.t()));
            //2p2h
            v_pphh.at(lmd2) += mult->dgemm(int_Cpppp_Vpphh, C_hhhh.t());
            v_pphh.at(lmd2) += mult->dgemm(int_Cpphh_Vpphht, C_hhpp.t());
            //            new_pphh.at(lmd2) += int_Cppph_Vphph * C_hhph.t();
            //            new_pphh.at(lmd2) -= int_Cppph_Vphph * C_hhhp.t();
            v_pphh.at(lmd2) += mult->dgemm(int_Cppph_Vphph, (C_hhph.t() - C_hhhp.t()));
            //            new_pphh.at(lmd2) -= int_Cpphp_Vphph * C_hhph.t();
            //            new_pphh.at(lmd2) += int_Cpphp_Vphph * C_hhhp.t();
            v_pphh.at(lmd2) += mult->dgemm(int_Cpphp_Vphph, (C_hhhp.t() - C_hhph.t()));
            //3p1h
            //            new_pphh.at(lmd2) += int_Cpppp_Vppph * C_hhph.t();
            //            new_pphh.at(lmd2) -= int_Cpppp_Vppph * C_hhhp.t();
            v_pphh.at(lmd2) += mult->dgemm(int_Cpppp_Vppph, (C_hhph.t() - C_hhhp.t()));
            v_pphh.at(lmd2) += mult->dgemm(int_Cppph_Vpppht, C_hhpp.t());
            v_pphh.at(lmd2) -= mult->dgemm(int_Cpphp_Vpppht, C_hhpp.t());
            //4p0h
            v_pphh.at(lmd2) += mult->dgemm(int_Cpppp_Vpppp, C_hhpp.t());
        }
        //phph
        {
            //0p4h
            v_phph.at(lmd2) = mult->dgemm(int_Cphhh_Vhhhh, C_phhh.t());
            //1p3h
            v_phph.at(lmd2) += mult->dgemm(int_Cphph_Vphhh, C_phhh.t());
            v_phph.at(lmd2) -= mult->dgemm(int_Cphhp_Vphhh, C_phhh.t());
            //            new_phph.at(lmd2) += int_Cphhh_Vphhht * C_phph.t();
            //            new_phph.at(lmd2) -= int_Cphhh_Vphhht * C_phhp.t();
            v_phph.at(lmd2) += mult->dgemm(int_Cphhh_Vphhht, (C_phph.t() - C_phhp.t()));
            //2p2h
            v_phph.at(lmd2) += mult->dgemm(int_Cphpp_Vpphh, C_phhh.t());
            v_phph.at(lmd2) += mult->dgemm(int_Cphhh_Vpphht, C_phpp.t());
            //            new_phph.at(lmd2) += int_Cphph_Vphph * C_phph.t();
            //            new_phph.at(lmd2) -= int_Cphph_Vphph * C_phhp.t();
            v_phph.at(lmd2) += mult->dgemm(int_Cphph_Vphph, (C_phph.t() - C_phhp.t()));
            //            new_phph.at(lmd2) -= int_Cphhp_Vphph * C_phph.t();
            //            new_phph.at(lmd2) += int_Cphhp_Vphph * C_phhp.t();
            v_phph.at(lmd2) += mult->dgemm(int_Cphhp_Vphph, (C_phhp.t() - C_phph.t()));
            //3p1h
            //            new_phph.at(lmd2) += int_Cphpp_Vppph * C_phph.t();
            //            new_phph.at(lmd2) -= int_Cphpp_Vppph * C_phhp.t();
            v_phph.at(lmd2) += mult->dgemm(int_Cphpp_Vppph, (C_phph.t() - C_phhp.t()));
            v_phph.at(lmd2) += mult->dgemm(int_Cphph_Vpppht, C_phpp.t());
            v_phph.at(lmd2) -= mult->dgemm(int_Cphhp_Vpppht, C_phpp.t());
            //4p0h
            v_phph.at(lmd2) += mult->dgemm(int_Cphpp_Vpppp, C_phpp.t());
        }
        //ppph
        {
            //0p4h
            v_ppph.at(lmd2) = mult->dgemm(int_Cpphh_Vhhhh, C_phhh.t());
            //1p3h
            v_ppph.at(lmd2) += mult->dgemm(int_Cppph_Vphhh, C_phhh.t());
            v_ppph.at(lmd2) -= mult->dgemm(int_Cpphp_Vphhh, C_phhh.t());
            //            new_ppph.at(lmd2) += int_Cpphh_Vphhht * C_phph.t();
            //            new_ppph.at(lmd2) -= int_Cpphh_Vphhht * C_phhp.t();
            v_ppph.at(lmd2) += mult->dgemm(int_Cpphh_Vphhht, (C_phph.t() - C_phhp.t()));
            //2p2h
            v_ppph.at(lmd2) += mult->dgemm(int_Cpppp_Vpphh, C_phhh.t());
            v_ppph.at(lmd2) += mult->dgemm(int_Cpphh_Vpphht, C_phpp.t());
            //            new_ppph.at(lmd2) += int_Cppph_Vphph * C_phph.t();
            //            new_ppph.at(lmd2) -= int_Cppph_Vphph * C_phhp.t();
            v_ppph.at(lmd2) += mult->dgemm(int_Cppph_Vphph, (C_phph.t() - C_phhp.t()));
            //            new_ppph.at(lmd2) -= int_Cpphp_Vphph * C_phph.t();
            //            new_ppph.at(lmd2) += int_Cpphp_Vphph * C_phhp.t();
            v_ppph.at(lmd2) += mult->dgemm(int_Cpphp_Vphph, (C_phhp.t() - C_phph.t()));
            //3p1h
            //            new_ppph.at(lmd2) += int_Cpppp_Vppph * C_phph.t();
            //            new_ppph.at(lmd2) -= int_Cpppp_Vppph * C_phhp.t();
            v_ppph.at(lmd2) += mult->dgemm(int_Cpppp_Vppph, (C_phph.t() - C_phhp.t()));
            v_ppph.at(lmd2) += mult->dgemm(int_Cppph_Vpppht, C_phpp.t());
            v_ppph.at(lmd2) -= mult->dgemm(int_Cpphp_Vpppht, C_phpp.t());
            //4p0h
            v_ppph.at(lmd2) += mult->dgemm(int_Cpppp_Vpppp, C_phpp.t());
        }
        //pppp
        {
            //0p4h
            v_pppp.at(lmd2) = mult->dgemm(int_Cpphh_Vhhhh, C_pphh.t());
            //1p3h
            v_pppp.at(lmd2) += mult->dgemm(int_Cppph_Vphhh, C_pphh.t());
            v_pppp.at(lmd2) -= mult->dgemm(int_Cpphp_Vphhh, C_pphh.t());
            //            new_pppp.at(lmd2) += int_Cpphh_Vphhht * C_ppph.t();
            //            new_pppp.at(lmd2) -= int_Cpphh_Vphhht * C_pphp.t();
            v_pppp.at(lmd2) += mult->dgemm(int_Cpphh_Vphhht, (C_ppph.t() - C_pphp.t()));
            //2p2h
            v_pppp.at(lmd2) += mult->dgemm(int_Cpppp_Vpphh, C_pphh.t());
            v_pppp.at(lmd2) += mult->dgemm(int_Cpphh_Vpphht, C_pppp.t());
            //            new_pppp.at(lmd2) += int_Cppph_Vphph * C_ppph.t();
            //            new_pppp.at(lmd2) -= int_Cppph_Vphph * C_pphp.t();
            v_pppp.at(lmd2) += mult->dgemm(int_Cppph_Vphph, (C_ppph.t() - C_pphp.t()));
            //            new_pppp.at(lmd2) -= int_Cpphp_Vphph * C_ppph.t();
            //            new_pppp.at(lmd2) += int_Cpphp_Vphph * C_pphp.t();
            v_pppp.at(lmd2) += mult->dgemm(int_Cpphp_Vphph, (C_pphp.t() - C_ppph.t()));
            //3p1h
            //            new_pppp.at(lmd2) += int_Cpppp_Vppph * C_ppph.t();
            //            new_pppp.at(lmd2) -= int_Cpppp_Vppph * C_pphp.t();
            v_pppp.at(lmd2) += mult->dgemm(int_Cpppp_Vppph, (C_ppph.t() - C_pphp.t()));
            v_pppp.at(lmd2) += mult->dgemm(int_Cppph_Vpppht, C_pppp.t());
            v_pppp.at(lmd2) -= mult->dgemm(int_Cpphp_Vpppht, C_pppp.t());
            //4p0h
            v_pppp.at(lmd2) += mult->dgemm(int_Cpppp_Vpppp, C_pppp.t());
        }
    }


    //Find f_hh
    for (int i = 0; i < nH; i++)
        for (int j = 0; j < nH; j++)
        {
            double u_ij = 0;
            for (int k = 0; k < nH; k++)
                u_ij += v_elem(i, k, j, k);
            f_hh(i, j) = h0(i, j) + u_ij;
        }
    //Find f_ph 
    for (int a = 0; a < nP; a++)
        for (int j = 0; j < nH; j++)
        {
            double u_aj = 0;
            for (int k = 0; k < nH; k++)
                u_aj += v_elem(a + nH, k, j, k);
            f_ph(a, j) = h0(a + nH, j) + u_aj;
        }
    //Find f_pp
    for (int a = 0; a < nP; a++)
        for (int b = 0; b < nP; b++)
        {
            double u_ab = 0;
            for (int k = 0; k < nH; k++)
                u_ab += v_elem(a + nH, k, b + nH, k);
            f_pp(a, b) = h0(a + nH, b + nH) + u_ab;
        }

    cout << "Transformation of Basis took " << timer.toc() << "s.\n";

    return;
}


