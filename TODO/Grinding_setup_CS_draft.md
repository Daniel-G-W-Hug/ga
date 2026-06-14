Let's refine the plan a little bit: The overarching target is to reproduce the main
results of the paper Tao et.al.. Basically, we want to reproduce the main results, i.e.
reproduce some of the diagrams shown in the paper. We can do that and visualize it using
ga_view (once we have created the results in the test files w/o graphics) by projecting
some of the 3DP results in the base plane which will be e423_3d (plane in which the chuck
and thus the wafer rotates) or in the plane e431_3dp of the global system defined below
(which shows the relative orientation of wafer and grinding tool).

In order to implement that we first need to reproduce the geometry and the kinematic
system as close as possible. That implies that we reproduce our results based on the
coordinate systems shown in fig. 1 of the paper: The root system is located on top of the
axis of the chuck (which holds the wafer on top). The the root system (=chuck surface)
is stationary with x (e1) pointing to the top (normal to the plane of the chuck that holds
the wafer later on), z (e3) is pointing to the right, and y (e2) is pointing to the front
(i.e. to the observer) to form a right-handed system. It is our intertial reference
system (-> "chuck_ctr_stat")

Then we basically have two independent kinematic chains build on top of this system (or
more if we want to introduce other sensor signals for example). All will be relative to
the root system.

The first chain start with a system that coincides with "chuck_ctr_stat" at t=0, but is
rotating with the chuck, so with rpm n_s (typ. 200-300 rpm) relative to the X axis of the
root system (-> "chuck_ctr_rot").

The second system of the first chain is shifted in X direction by the average wafer
thickness (tw_avg) to the positive X direcation. Its yz-plane represents the average wafer
thickness. Any deviation from X == 0 represents the current surface profile of the top
side of the wafer vs. the average. So it is located at the average top surface of the
wafer. This might change over time if tw_avg changes, i.e. when the wafer becomes thinner,
but currently it is not clear on which time scale that will happen. At beginning it can be
regarded nearly constant (i.e. does not change significantly before the chuck itself has
seen 10 revolutions, intermediate changes are covered in the profile). The three
directions coincide with the directions of the global system at t=0 and the Y and Z
directions rotate with its parent system, i.e. all directions continuously coincide with
the parent system (-> "wafer_top_avg_rot").

The second chain is also rooted at "chuck_ctr_stationary", but now describes the rotating
spindle and the grinding wheel attached to it. Later we will be interested in relative
speed at the outer diameter of the grinding wheel. The relative speed at the position
where the wafer surface meets the surface of the grinding tool will be an important input
to the grinding model.

The first system of the second chain is the center of mass of the spindle. Its x-direction
coincides with the z-direction of the root system, and the y-direction coincides with the
y-direction of the root system. The z-direction is directed in the negative x-direction of
the root system (-> "spindle_cm_stat"). The spindle pose has all degrees of freedom vs.
the root system. Axial motion is along its z-direction, radial motion in x- and/or
y-direction i.e. a translational shift of the rotation axis the spindle vs. the rotation
axis of the chuck. A tilting motion is a rotational shift of the rotation axis of the
spindle vs. the rotation axis of the chuck.

Attached to the spindle we have the grinding tool, which has a displacement in positive
z-direction of the spindle system. It represents the distance of the average surface of
the grinding tool from the center of gravity of the spindle (l3 in fig. 1). This
displacement might change over time, if the tool shows some wear. To begin with we can
assume that it only changes slowly, i.e. with small temporal dynamic. The system coincides
at t=0 with "spindle_cm_stat", but rotates around the z-axis with rpm n_s (typ. 3000 rpm).
The system shows specific positions z = l1 and z=l2 where damping and spring elements will
act, but that can be covered once the dynamics is modelled. For now only the distance to
the center of mass is relevant (-> "tool_top_avg_rot").

A third system in that second chain will be located at the outer diameter of the grinding
wheel. It is located at the position x representing the diameter of the wheel and it
rotates with the grinding wheel, i.e. all axis coincide with the parent directions, but
the whole system is shifted by the radius of the grinding wheel to the x-direction of the
parent system (-> "tool_surface_avg_at_R").

Does that description match with your interpretation of fig 1? If yes, please add that to
the plan in order to set up the kinematic system accordingly.

One early try to use the system and show it to be functional should be to reproduce the
curves shown a gray dashed lines in fig. 1 or as red dashed lines in fig. 7.They show the
trajectory of a grain in the tool surface shown in the system of the rotating wafer. Those
patterns can be shown for specific combinations of n_s, n_w and the corresponding wafer
and tool diameters. That can be done before we go into the dynamics. ga_view can be used
for visualization.

It would be good when the initial trials we do to replicate the Sommerfeld effect already
use a similar setup, so that we can get used to the complex system step by step.