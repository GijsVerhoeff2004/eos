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
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

using namespace std::literals::complex_literals;

namespace eos
{

    template <> struct Implementation<BToDPseudoscalar>
    {
            std::shared_ptr<Model>                        model;
            QuarkFlavorOption                             opt_q;
            QuarkFlavorOption                             opt_D;
            LightMesonOption                              opt_p;
            UsedParameter                                 hbar;
            UsedParameter                                 tau;
            UsedParameter                                 mB;
            UsedParameter                                 mD;
            UsedParameter                                 mP;
            RestrictedOption                              opt_rep;
            std::shared_ptr<NonleptonicAmplitudes<PToDP>> nl_amplitudes;
            std::shared_ptr<NonleptonicAmplitudes<PToDP>> cp_nl_amplitudes;
            std::shared_ptr<NonleptonicAmplitudes<PToDP>> Bbar_nl_amplitudes;

            // B mixing
            std::shared_ptr<UsedParameter>                life_time_difference;
            std::function<double()>                       yq;
            std::function<double()>                       phiB;
            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"_ov), p, o)),
            opt_q(o, options, "q"_ok),
            opt_D(o, options, "D"_ok),
            opt_p(o, options, "P"_ok),
            hbar(p["QM::hbar"], u),
            tau(p["life_time::B_" + opt_q.str()], u),
            mB(p["mass::B_" + opt_q.str()], u),
            mD(p["mass::D_" + opt_D.str()], u),
            mP(p["mass::" + opt_p.str()], u),
            opt_rep(o, options, "representation"_ok),
            nl_amplitudes(NonleptonicAmplitudeFactory<PToDP>::create("B->DP::" + opt_rep.value(), p, o + Options{{"cp-conjugate"_ok, "false"_ov}})),
            cp_nl_amplitudes(NonleptonicAmplitudeFactory<PToDP>::create("B->DP::" + opt_rep.value(), p, o + Options{{"cp-conjugate"_ok, "true"_ov}})),
            Bbar_nl_amplitudes(NonleptonicAmplitudeFactory<PToDP>::create("B->DP::" + opt_rep.value(), p, o + Options{{"cp-conjugate"_ok, "false"_ov}} + Options{{"B-bar"_ok, "true"_ov}}))
            {
                Context ctx("When constructing B->DP observable");

                switch (opt_q.value())
                {
                    case QuarkFlavor::up:
                        life_time_difference = nullptr;
                        yq                   = []() { return 0.0; };
                        phiB                 = []() { return 0.0; };
                        break;

                    case QuarkFlavor::down:
                        life_time_difference = std::make_shared<UsedParameter>(p["life_time::Delta_B_d"], u);
                        yq                   = [this]() { return this->life_time_difference->evaluate() / (2.0 / tau); };
                        phiB                 = [this]() { return 2.0 * arg(model->ckm_tb() * conj(model->ckm_td())); };
                        break;

                    case QuarkFlavor::strange:
                        life_time_difference = std::make_shared<UsedParameter>(p["life_time::Delta_B_s"], u);
                        yq                   = [this]() { return this->life_time_difference->evaluate() / (2.0 / tau); };
                        phiB                 = [this]() { return 2.0 * arg(model->ckm_tb() * conj(model->ckm_ts())); };
                        break;

                    default: throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
                }

                u.uses(*model);
                u.uses(*nl_amplitudes);
                u.uses(*cp_nl_amplitudes);
                u.uses(*Bbar_nl_amplitudes);
            }

            double
            decay_width() const
            {
                const double mB  = this->mB();
                const double mD  = this->mD();
                const double mP  = this->mP();
                const double pre = 1.0 / 16.0 / M_PI / power_of<3>(mB) * std::sqrt(lambda(mB * mB, mD * mD, mP * mP));

                return pre * std::norm(nl_amplitudes->amplitude());
            }

            double
            cp_decay_width() const
            {
                const double mB  = this->mB();
                const double mD  = this->mD();
                const double mP  = this->mP();
                const double pre = 1.0 / 16.0 / M_PI / power_of<3>(mB) * std::sqrt(lambda(mB * mB, mD * mD, mP * mP));

                return pre * std::norm(cp_nl_amplitudes->amplitude());
            }

            double
            mixing_induced_cp_asymmetry() const
            {
                const complex<double> amp      = nl_amplitudes->amplitude();
                const complex<double> Bbar_amp = Bbar_nl_amplitudes->amplitude();

                // Assumes the mixing parameter ratio q / p to be a pure phase
                const complex<double> xif = -std::exp(-1.0i * phiB()) * Bbar_amp / amp;

                return 2 * std::imag(xif) / (1 + std::norm(xif));
            }

            double
            a_Delta_Gamma() const
            {
                const complex<double> amp      = nl_amplitudes->amplitude();
                const complex<double> Bbar_amp = Bbar_nl_amplitudes->amplitude();

                // Assumes the mixing parameter ratio q / p to be a pure phase
                const complex<double> xif = -std::exp(-1.0i * phiB()) * Bbar_amp / amp;

                return 2 * std::real(xif) / (1 + std::norm(xif));
            }
    };

    const std::vector<OptionSpecification> Implementation<BToDPseudoscalar>::options{
        Model::option_specification(),
        NonleptonicAmplitudeFactory<PToDP>::option_specification(),
        {              "q"_ok,                                                                                            { "u"_ov, "d"_ov, "s"_ov } },
        {              "D"_ok,                                                                                            { "u"_ov, "d"_ov, "s"_ov } },
        {              "P"_ok, { "pi^0"_ov, "pi^+"_ov, "pi^-"_ov, "K_d"_ov, "Kbar_d"_ov, "K_S"_ov, "K_u"_ov, "Kbar_u"_ov, "eta"_ov, "eta_prime"_ov } },
        { "representation"_ok,                                                                                                         { "QCDF"_ov } }
    };

    BToDPseudoscalar::BToDPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDPseudoscalar>(new Implementation<BToDPseudoscalar>(parameters, options, *this))
    {
    }

    BToDPseudoscalar::~BToDPseudoscalar() {}

    double
    BToDPseudoscalar::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BToDPseudoscalar::branching_ratio() const
    {
        return _imp->decay_width() * _imp->tau() / _imp->hbar();
    }

    double
    BToDPseudoscalar::cp_branching_ratio() const
    {
        return _imp->cp_decay_width() * _imp->tau() / _imp->hbar();
    }

    double
    BToDPseudoscalar::avg_branching_ratio() const
    {
        return 0.5 * (branching_ratio() + cp_branching_ratio());
    }

    double
    BToDPseudoscalar::exp_branching_ratio() const
    {
        const double yq = _imp->yq();
        return avg_branching_ratio() * (1 + _imp->a_Delta_Gamma() * yq) / (1 - power_of<2>(yq));
    }

    double
    BToDPseudoscalar::cp_asymmetry() const
    {
        return (branching_ratio() - cp_branching_ratio()) / (branching_ratio() + cp_branching_ratio());
    }

    double
    BToDPseudoscalar::mixing_induced_cp_asymmetry() const
    {
        return _imp->mixing_induced_cp_asymmetry();
    }

    double
    BToDPseudoscalar::a_Delta_Gamma() const
    {
        return _imp->a_Delta_Gamma();
    }

    const std::string BToDPseudoscalar::description = "\
    The decay B->DP, where D is a charmed and P a light pseudoscalar. The option 'D' selects the light antiquark flavour of the charmed meson.";

    const std::set<ReferenceName> BToDPseudoscalar::references{
        "BBNS:2000A"_rn,
    };

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::begin_options()
    {
        return Implementation<BToDPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::end_options()
    {
        return Implementation<BToDPseudoscalar>::options.cend();
    }
} // namespace eos
