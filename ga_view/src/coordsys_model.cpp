#include "coordsys_model.hpp"

[[maybe_unused]] int Coordsys_model::add_pt(pt2d const& pt_in, pt2d_mark const m)
{

    pt.push_back(pt_in);
    pt_mark.push_back(m);

    mark_id new_id;
    new_id.id = unique_id++;
    pt_id.push_back(new_id);

    return new_id.id;
}

//
// hint: using ln2d = std::vector<pt2d>;
//
[[maybe_unused]] int Coordsys_model::add_ln(ln2d const& vp_in, ln2d_mark const m)
{

    // the separate copy should not be needed, since it is done in push_back
    // anyway
    //
    // std::vector<pt2d> v;
    // std::copy(vp_in.begin(), vp_in.end(), std::back_inserter(v));
    //
    ln.push_back(vp_in);
    ln_mark.push_back(m);

    mark_id new_id;
    new_id.id = unique_id++;
    ln_id.push_back(new_id);

    if (m.mark_pts == true) { // add points of line to pts marked in model

        for (int i = 0; i < vp_in.size(); i += m.delta) {

            pt.push_back(vp_in[i]);
            pt_mark.push_back(m.pm);

            mark_id new_pt_id;
            new_pt_id.id = unique_id++;
            new_pt_id.linked_to_id = new_id.id;
            pt_id.push_back(new_pt_id);
        }
    }

    return new_id.id;
}

[[maybe_unused]] int Coordsys_model::add_vt(vt2d const& vt_in, vt2d_mark const m)
{

    vt.push_back(vt_in);
    vt_mark.push_back(m);

    mark_id new_id;
    new_id.id = unique_id++;
    vt_id.push_back(new_id);

    return new_id.id;
}

[[maybe_unused]] int Coordsys_model::add_apt(pt2d const& pt_in)
{

    apt.push_back(pt_in);

    mark_id new_id;
    new_id.id = unique_id++;
    apt_id.push_back(new_id);

    return apt.size() - 1;
}

[[maybe_unused]] int Coordsys_model::add_avt(avt2d const& avt_in)
{

    avt.push_back(avt_in);

    mark_id new_id;
    new_id.id = unique_id++;
    avt_id.push_back(new_id);

    return avt.size() - 1;
}

[[maybe_unused]] int Coordsys_model::add_abivt(abivt2d const& abivt_in)
{

    abivt.push_back(abivt_in);

    mark_id new_id;
    new_id.id = unique_id++;
    avt_id.push_back(new_id);

    return abivt.size() - 1;
}

[[maybe_unused]] int Coordsys_model::add_aproj(aproj2d const& aproj_in)
{

    aproj.push_back(aproj_in);

    mark_id new_id;
    new_id.id = unique_id++;
    avt_id.push_back(new_id);

    return abivt.size() - 1;
}

[[maybe_unused]] int Coordsys_model::add_arefl(arefl2d const& arefl_in)
{

    arefl.push_back(arefl_in);

    mark_id new_id;
    new_id.id = unique_id++;
    avt_id.push_back(new_id);

    return abivt.size() - 1;
}

void Coordsys_model::clear()
{
    unique_id = 0;

    pt.clear();
    pt_mark.clear();
    pt_id.clear();

    ln.clear();
    ln_mark.clear();
    ln_id.clear();

    m_label.clear();
}