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
 * File:   HFbasis.h
 * Author: toffyrn
 *
 * Created on 16. februar 2012, 17:21
 */

#ifndef HFBASIS_H
#define	HFBASIS_H

#include "config.h"

#include "System.h"
#include "GEMM.h"

namespace toffyrn
{
    namespace libUNK
    {

        /**
         * @brief HFsys is a system with a transformed basis. All elements are rewritten
         * as using a unitary matrix C obtained by the Hartree Fock method.
         * 
         * The new basis is then
         * \f$ |\psi \rangle = \sum_{\gamma} C_{\psi, \gamma} |\gamma \rangle \f$.
         * In this sense all matrix elements get new values;
         * \f$ h^{(0)}_{pq} = \sum_{\alpha \beta} C_{p\alpha} C_{q\beta} h^{(0)}_{\alpha \beta} \f$
         * and 
         * \f$ \langle pq||rs \rangle = \sum_{\alpha \beta} \sum_{\gamma \delta}
         * C_{p\alpha} C_{q\beta} C_{r\gamma} C_{s\delta} \langle \alpha \beta || \gamma \delta \rangle \f$
         * 
         * The Basis object is copied, as the channels and configurations are 
         * unaltered.
         */
        class HFsys : public System
        {
        public:

            /**
             * @brief Create a new system similar to originalSys, but with transformed elements.
             * @param originalSys Original system to copy Basis and interactions from.
             * @param Coeff Coefficient matrix from HF, used to transform interactions.
             * @param mult Optionally, a matrix multiplicator. 
             */
            HFsys(System const * originalSys, arma::mat const &Coeff, GEMM * mult = NULL);

            /**
             * Destructor, deleting matrix multiplicator "mult".
             */
            virtual ~HFsys();

            virtual double f_elem(std::size_t p, std::size_t q) const;
            virtual double v_elem(
                    std::size_t p, std::size_t q,
                    std::size_t r, std::size_t s) const;

        protected:
            /**
             * @brief Transformation of elements is done in this method.
             * @param C coefficient matrix.
             * @param origSys Original System
             */
            virtual void transformElements(arma::mat const &C, System const * origSys);
        private:
            /** Matrix multiplicator */
            GEMM * mult;

            /** Standard matrix multiplicator */
            GEMM mult_standard;

        };

    }
}
#endif	/* HFBASIS_H */

