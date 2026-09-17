/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/b-decays/b-to-d-psd.hh>
#include <eos/maths/complex.hh>
#include <eos/observable.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class BToDPseudoscalarTest : public TestCase
{
    public:
        BToDPseudoscalarTest() :
            TestCase("b_to_d_psd_test")
        {
        }

        virtual void
        run() const
        {
            static const double eps = 1.0e-5;

            Parameters p                         = Parameters::Defaults();
            p["nonleptonic::Re{alpha1}@QCDF-DP"] = 1.0;
            p["nonleptonic::Re{alpha2}@QCDF-DP"] = 0.2;
            p["nonleptonic::Re{b1}@QCDF-DP"]     = 0.01;

            // colour-allowed emission with a single weak phase
            {
                Options o{
                    { "representation"_ok, "QCDF"_ov },
                    {              "q"_ok,    "d"_ov },
                    {              "D"_ok,    "d"_ov },
                    {              "P"_ok, "pi^-"_ov },
                };

                BToDPseudoscalar d(p, o);

                TEST_CHECK(d.branching_ratio() > 0.0);
                // all coefficients share the CKM factor V_cb V_ud^*, hence no direct CP asymmetry arises
                TEST_CHECK_NEARLY_EQUAL(d.cp_asymmetry(), 0.0, eps);
                TEST_CHECK_RELATIVE_ERROR(d.avg_branching_ratio(), d.branching_ratio(), eps);
            }

            // colour-suppressed emission
            {
                Options o{
                    { "representation"_ok, "QCDF"_ov },
                    {              "q"_ok,    "d"_ov },
                    {              "D"_ok,    "u"_ov },
                    {              "P"_ok, "pi^0"_ov },
                };

                BToDPseudoscalar d(p, o);

                TEST_CHECK(d.branching_ratio() > 0.0);
            }

            // charge is not conserved in Bbar^0 -> D^+ pi^0, hence the amplitude vanishes
            {
                Options o{
                    { "representation"_ok, "QCDF"_ov },
                    {              "q"_ok,    "d"_ov },
                    {              "D"_ok,    "d"_ov },
                    {              "P"_ok, "pi^0"_ov },
                };

                BToDPseudoscalar d(p, o);

                TEST_CHECK_NEARLY_EQUAL(d.decay_width(), 0.0, eps);
            }

            // vanishing QCDF coefficients yield a vanishing decay width
            {
                Parameters q                         = Parameters::Defaults();
                q["nonleptonic::Re{alpha1}@QCDF-DP"] = 0.0;

                Options o{
                    { "representation"_ok, "QCDF"_ov },
                    {              "q"_ok,    "d"_ov },
                    {              "D"_ok,    "d"_ov },
                    {              "P"_ok, "pi^-"_ov },
                };

                BToDPseudoscalar d(q, o);

                TEST_CHECK_NEARLY_EQUAL(d.decay_width(), 0.0, eps);
            }

            // the registered observables agree with the class interface
            {
                Options o{
                    { "representation"_ok, "QCDF"_ov },
                };

                ObservablePtr obs = Observable::make("B^0->D^+pi^-::BR", p, Kinematics{}, o);

                TEST_CHECK(ObservablePtr() != obs);

                Options oc{
                    { "representation"_ok, "QCDF"_ov },
                    {              "q"_ok,    "d"_ov },
                    {              "D"_ok,    "d"_ov },
                    {              "P"_ok, "pi^-"_ov },
                };

                BToDPseudoscalar d(p, oc);

                TEST_CHECK_RELATIVE_ERROR(obs->evaluate(), d.branching_ratio(), eps);
            }
        }
} b_to_d_psd_test;
