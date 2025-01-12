//
// author: Daniel Hug, 2024
//

#include "coordsys_model.hpp"

[[maybe_unused]] size_t Coordsys_model::add_pt(pt2d const& pt_in, pt2d_mark const& m)
{

    pt.push_back(pt_in);
    pt_mark.push_back(m);

    return pt.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_pt(pt2dp const& pte_in, pt2d_mark const& m)
{

    pte.push_back(pte_in);
    pte_mark.push_back(m);

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
// hint: using ln2de = std::vector<pt2dp>;
//
[[maybe_unused]] size_t Coordsys_model::add_ln(ln2de const& vpe_in, ln2d_mark const& m)
{
    lne.push_back(vpe_in);
    lne_mark.push_back(m);

    if (m.mark_pts == true) { // add points of line to pts marked in model

        for (size_t i = 0; i < vpe_in.size(); i += m.delta) {

            pte.push_back(vpe_in[i]);
            pte_mark.push_back(m.pm);
        }
    }

    return lne.size() - 1;
}

[[maybe_unused]] size_t Coordsys_model::add_vt(vt2d const& vt_in, vt2d_mark const& m)
{

    vt.push_back(vt_in);
    vt_mark.push_back(m);

    return vt.size() - 1;
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

void Coordsys_model::clear()
{
    pt.clear();
    pt_mark.clear();

    pte.clear();
    pte_mark.clear();

    ln.clear();
    ln_mark.clear();

    lne.clear();
    lne_mark.clear();

    vt.clear();
    vt_mark.clear();

    apt.clear();
    avt.clear();
    abivt.clear();
    aproj.clear();
    arefl.clear();

    m_label.clear();
}