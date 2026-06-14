// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "coordsys_model.hpp"

[[maybe_unused]] size_t Coordsys_model::add_pt(pt2d const& pt_in, pt2d_mark const& m)
{

    pt.push_back(pt_in);
    pt_mark.push_back(m);

    return pt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_pt(pt2dp const& ptp_in, pt2d_mark const& m)
{

    ptp.push_back(ptp_in);
    ptp_mark.push_back(m);

    return pt.size() - 1;
}

//
// hint: using ln2d = std::vector<pt2d>;
//
[[maybe_unused]] size_t Coordsys_model::add_ln(ln2d const& vp_in, ln2d_mark const& m)
{
    ln.push_back(vp_in);
    ln_mark.push_back(m);

    if (m.mark_pts == true) { // add points of line to pts marked in model

        for (size_t i = 0; i < vp_in.size(); i += m.delta) {

            pt.push_back(vp_in[i]);
            pt_mark.push_back(m.pm);
        }
    }

    return ln.size() - 1;
}

//
// hint: using cln2dp = std::vector<pt2dp>;
//
[[maybe_unused]] size_t Coordsys_model::add_ln(cln2dp const& vp_in, ln2d_mark const& m)
{
    clnp.push_back(vp_in);
    clnp_mark.push_back(m);

    if (m.mark_pts == true) { // add points of line to pts marked in model

        for (size_t i = 0; i < vp_in.size(); i += m.delta) {

            ptp.push_back(vp_in[i]);
            ptp_mark.push_back(m.pm);
        }
    }

    return clnp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_vt(vt2d const& vt_in, vt2d_mark const& m)
{

    vt.push_back(vt_in);
    vt_mark.push_back(m);

    return vt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_bivtp(bivt2dp const& bivtp_in,
                                                  bivt2dp_mark const& m)
{

    bivtp.push_back(bivtp_in);
    bivtp_mark.push_back(m);

    return bivtp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_apt(pt2d const& pt_in)
{

    apt.push_back(pt_in);

    return apt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_avt(avt2d const& avt_in)
{

    avt.push_back(avt_in);

    return avt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_abivt(abivt2d const& abivt_in)
{

    abivt.push_back(abivt_in);

    return abivt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_aproj(aproj2d const& aproj_in)
{

    aproj.push_back(aproj_in);

    return aproj.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_arefl(arefl2d const& arefl_in)
{

    arefl.push_back(arefl_in);

    return arefl.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_abivtp(abivt2dp const& abivtp_in)
{

    abivtp.push_back(abivtp_in);

    return abivtp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_areflp(arefl2dp const& areflp_in)
{

    areflp.push_back(areflp_in);

    return areflp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_akinem(akinematic2dp const& akinem_in)
{

    akinemp.push_back(akinem_in);

    return akinemp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_aode(aode_spring2d const& aode_in)
{

    aode.push_back(aode_in);

    return aode.size() - 1;
}

[[maybe_unused]] size_t
Coordsys_model::add_aode_plate(aode_plate_pga2dp const& aode_plate_in)
{
    aode_plate.push_back(aode_plate_in);

    return aode_plate.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_aframe_trafo(aframe_trafo const& aft_in)
{
    aft.push_back(aft_in);

    return aft.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_merry_go_round(amerry_go_round const& amgr_in)
{
    amgr.push_back(amgr_in);

    return amgr.size() - 1;
}

[[maybe_unused]] size_t
Coordsys_model::add_double_pendulum(adouble_pendulum const& adp_in)
{
    adp.push_back(adp_in);

    return adp.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_four_bar(afour_bar const& afb_in)
{
    afb.push_back(afb_in);

    return afb.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_open_vs_closed(aopen_vs_closed const& aovc_in)
{
    aovc.push_back(aovc_in);

    return aovc.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_planar_delta(aplanar_delta const& apld_in)
{
    apld.push_back(apld_in);

    return apld.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_grinding_marks(agrinding_marks const& agm_in)
{
    agm.push_back(agm_in);

    return agm.size() - 1;
}

void Coordsys_model::clear()
{
    pt.clear();
    pt_mark.clear();

    ptp.clear();
    ptp_mark.clear();

    ln.clear();
    ln_mark.clear();

    clnp.clear();
    clnp_mark.clear();

    vt.clear();
    vt_mark.clear();

    bivtp.clear();
    bivtp_mark.clear();

    apt.clear();
    avt.clear();
    abivt.clear();
    aproj.clear();
    arefl.clear();
    abivtp.clear();
    areflp.clear();
    akinemp.clear();
    aode.clear();
    aode_plate.clear();
    amgr.clear();
    adp.clear();
    afb.clear();
    aovc.clear();
    apld.clear();
    agm.clear();
    aft.clear();

    m_label.clear();
}