#include "coordsys_model.hpp"

[[maybe_unused]] int Coordsys_model::add_pt(const pt2d& p_in, const pt2d_mark m)
{

    pt.push_back(p_in);
    pt_mark.push_back(m);

    mark_id new_id;
    new_id.id = unique_id++;
    pt_id.push_back(new_id);

    return new_id.id;
}

//
// hint: using ln2d = std::vector<pt2d>;
//
[[maybe_unused]] int Coordsys_model::add_ln(const ln2d& vp_in, const ln2d_mark m)
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

[[maybe_unused]] int Coordsys_model::add_vt(const vt2d& v_in, const vt2d_mark m)
{

    vt.push_back(v_in);
    vt_mark.push_back(m);

    mark_id new_id;
    new_id.id = unique_id++;
    vt_id.push_back(new_id);

    return new_id.id;
}

[[maybe_unused]] int Coordsys_model::add_apt(const pt2d& p_in)
{

    apt.push_back(p_in);

    mark_id new_id;
    new_id.id = unique_id++;
    apt_id.push_back(new_id);

    return new_id.id;
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