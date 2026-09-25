# Changes

- 2024/02: add first routines for geometric algebra (GA)
- 2024/04: complete geometric products for 2d and 3d EGA incl. ga_test
- 2024/05: add ga_lua for scripting capabilities
- 2024/07: add ga_view as qt-based viewer for 2d applications with interactive point
           movements and consistent graphicsview system for active/passive items
- 2024/08: moved EGA (Euclidean GA) in own namespace and created PGA (Projective GA);
           implemented comprehensive type system redesign with strong typing (scalar_t);
           added support for up to 4D EGA algebras; added Windows/MSVC compatibility;
           optimized sandwich product rotations; implemented full duality support
           (Macdonald default, Hestenes optional)
- 2024/12: reworked structure of ega2d and ega3d:
                - common base types established
                - missing products added: dot, wdg, rwdg, inner, fat_dot,
                  left- and right complements
                - duality and complements added
                - unified nrm functions for all types
- 2025/01: added ga_prdxpr to automatically generate product expressions for all products;
           implemented comprehensive PGA operations including contractions, bulk/weight
           functionality, and projective transformations; added inverse operations
           (gr_inv, rev, conj) for PGA; enhanced robustness of active elements in ga_view
- 2025/02: started to implement parser for string manipulation with ga_prdxpr to
           manipulate and simplify expressions resulting from sandwich products in order
           to optimize performance (still needs completion for automatic simplification)
           Reworked pga2dp and implemented pga3dp based on expressions generated using
           ga_prdxpr
- 2025/03: implemented pga2dp and pga3dp including tests; reworked everything into  
           definitions as proposed by Lengyel consistently for ega2d, ega3d, pga2dp,
           pga3dp; added regressive geometric products (rgpr) and motors for PGA;
           implemented bulk/weight dual operations; enhanced expression parsing and
           tree visualization in ga_prdxpr; started to implement latex documentation
- 2025/05: implemented operations and visualization of rotors in ega and motors in pga;
           added commutator products (asymmetric part of geometric product); enhanced
           reflections and projective transformations; corrected half-angle formulas
           for PGA cases; improved active reflection documentation and visualization
- 2025/06: implemented unified error handling for division by zero; optimized split
           between files in ga library for improved modularity; reworked ga_prdxpr for
           more modularized math expression generation and preparation for automated
           transformation of sandwich product expressions to show they are fully
           equivalent to matrix * vector expression; improved modularization of ga_lua and
           ga_test with own files for each respective algebra
- 2025/07: enhanced ga_lua with interactive shell mode and readline support for better
           user experience; added comprehensive PGA (Projective Geometric Algebra) support
           to ga_lua including 2dp and 3dp operations; improved file structure with
           separate test files for each algebra; significant refactoring of ga_prdxpr
           with new generator-based implementation and automated transformation system
           for sandwich products
- 2025/08: implemented automatic generation of algebra rules system enabling easier
           extension with new algebras; added comprehensive congruence tests for all
           supported algebras (ega2d, ega3d, pga2dp, pga3dp); enhanced ga_view with
           open kinematic chain visualization and improved active reflections; unified
           color schemes for rotations; extensive code quality improvements including
           removal of compiler warnings, consistent zero checking, and improved type
           safety with disallowed implicit scalar conversions; added force visualization
           in 2D ga_view; implemented cross-platform modular dependency management system
           supporting flexible use of system packages vs FetchContent for better
           developer experience across macOS, Windows, and Linux
- 2025/09: enhanced documentation with metric signatures, bulk/weight concepts, and norms;
           added reverse commutator product (rcmt) consistently across EGA and PGA
           algebras in ga_prdxpr; implemented transwedge products for all algebras with
           application tests; introduced an application-specific overlap example;
           improved ga_lua with readline support, command-line option handling, and copy
           semantics for user-defined types; enhanced ga_view visualization with
           bivector/vector labels and force vector display; improved story-telling flow
           through reordered examples in ga_view; refined Windows/vcpkg integration;
           optimized fmt library integration using nested_formatter for proper GA type
           formatting; automatically generate algebra rules in ga_prdxpr
- 2025/10: restructured multivector type system eliminating cross-tag comparison
           requirements on GCC within the algebra-specific headers (mvec2d, mvec3d,
           mvec2dp, mvec3dp, mvec4d); reorganized test suite converting implementation
           files to header-only structure for better modularity; revised projection
           operations in pga2dp and pga3dp; enhanced documentation covering user-defined
           types, projective vs plane-based geometric algebra distinctions, and expanded
           literature database; added comprehensive demo scripts to ga_lua (contraction,
           reflection, rotation, inversion, expansion examples); further developed
           an application-specific example with improved overlap calculations and
           enhanced utility functions; test case added to show how to recover the perpendicular
           distance from a bivector created with the wedge product by applying the
           geometric product
- 2025/11: added dynamic systems support with spring-mass system implementation in
           pga2dp; extended motor operations for robotics applications with kinematic
           chains; enhanced physics modelling documentation including forces as lines
           (forque); improved ga_prdxpr with automatic dual calculation via extended
           metric and STA4D product tables; unified function naming (support(),
           is_congruent()); added DualNum2dp/3dp infrastructure with tests and lua
           integration; reworked transformations in pga2dp/pga3dp including screw motion;
           improved ga_view with connecting arrows for reflections and updated angle
           calculation; added bulk_normalize to lua; added metrics and regressive metrics
           to user constants; MSVC compatibility improvements
- 2025/12: added regressive wedge and geometric products for PGA multivectors in ga lib
           and ga_lua; implemented multivector normalization; added motion modelling
           chapter for EGA and PGA with tests in ega3d and pga3dp; extended commutator
           products for trivector-vector combinations; added conditional rotor/motor
           normalization; implemented product table output from ga_prdxpr for easy
           reference; enhanced documentation with complements, duals, dot products for
           pga3dp; added momentum tests and documentation; stricter numeric_types type
           checking; simplified rk_integration file structure; extended wdg/rwdg for
           even/odd PGA multivectors; fixed MSVC compatibility issues
- 2026/01: switched pga2dp basis bivectors from (e23, e31, e12) to (e31, e32, e12) for
           better compatibility with pga3dp where lines encode direction in bivectors;
           added metric signature as user input for each algebra in ga_prdxpr removing
           hard-coded elements; implemented recursive extended metric extraction from
           wedge/dot product rules with validation tests; added consistency checks for
           user input of multivector basis; updated rule generator to use header file
           input instead of hard-coding; added dynamics foundation with inertia map and
           bivector rate of change documentation; added comparison chapter of different
           approaches with table of contents; added pga2dp tables and basis variant
           comparison; renamed support() to sup(); added LICENSE file; prepared ideas
           for direct C++ codegen by ga_prdxpr; adapted ga_view for active bivectors
           after basis change
- 2026/02: introduced initial physics infrastructure for PGA: inertia map, inertia inverse,
           Runge-Kutta integration (rk_integration) with EGA3D reference test cases;
           added exp() and sqrt() for PGA with documentation and tests for the translation
           case; added get_motor_from_lines() for PGA 3D; added rgr() (right geometric
           reverse product) for PGA; refined PGA motor unitization; corrected reflection
           formulas for PGA; bulk_nrm_sq() and weight_nrm_sq() now return plain scalar
           values; added constexpr constructors and lv()/lm() helpers to bivec3dp;
           reordered movexdp() and added support for moving a motor in pga2dp;
           added bulk() and weight() to ga_lua for PGA; updated pga2dp output to new basis
           vectors; silenced unused-variable compiler warnings
- 2026/03: added exp() and sqrt() with corresponding tests for EGA; reworked and verified
           optimized implementations of rotate(), move2dp() and move3dp(); extended physics
           modelling with inertia for continuous mass distributions (flat plate 2D, cuboid
           3D) including documentation; extended coupled Newton-Euler equations and
           constrained motion documentation; implemented physics simulation for rigid body
           case with verification tests; added coupling between GA dynamics and classical
           dynamics equations in test cases; added combined rigid motion (translation and
           rotation) for a plate in 2D; implemented ODE-based plate simulation
           (active_plate_ode) in ga_view; added kinematics tests for 2D case; implemented
           pga3dp physics test for cuboid; added active frame transformation visualization
           to ga_view demonstrating body-frame (B_b) and world-frame (B_w) motion;
           added diagram_legend to ga_view; introduced separate executables for 2dp and 3dp
           application test cases; reordered ga_view test cases for better narrative flow;
           corrected ga_algebra.hpp for pga2dp basis vector change; removed shadowing
           warnings identified by cppcheck static analysis; switched selected tests from
           absolute to relative accuracy; added VTK 3D visualization architecture analysis
           and prep work; added GitHub repository setup guide for branch protection and
           coworker workflow
- 2026/04: implemented an initial complete Python wrapper for the ga library in two
           pieces: ga_bindgen (a libclang-based scanner that produces a manifest of the
           public C++ API) and ga_py (a nanobind extension generated from the manifest,
           exposing the full ga/ surface as ga_py.ega and ga_py.pga submodules); ported
           the Python wrapper to Windows; added vec3dp/vec2dp constructors from
           vec3d/vec2d; completed testing for physics modelling of rotations in moving
           systems; resolved bugs in several templated product expressions that surfaced
           during the wrapper work; established a Python test suite with 411 tests
           covering constants, grade lookup, EGA/PGA algebraic identities (randomized via
           hypothesis), and a JSON cross-check against the C++ reference output to ensure
           C++ and Python give the same results; renamed left_/right_ to l_/r_
           consistently across complements and duals for a uniform naming pattern; moved
           norms out of the fundamental type definitions into the algebras themselves;
           added regressive multivector products with ps and s for PGA; completed
           documentation on rotational accelerations; added active_ode_plate vector
           addition in ga_view; restructured ga_py / ga_bindgen documentation for
           first-time users (six-section ga_py/README with quick-start, demo script at
           ga_py/demo/hello_ga.py, with earlier scratch notes folded into the
           READMEs); various documentation improvements and a
           spellchecking pass
- 2026/07: added STA4D (space-time algebra, G(1,3,0)) with geometric-product rotor
           sandwich transforms, boosts, exp/log/sqrt, and closed-form transform_opt
           helpers; added the regressive inverse rinv() for PGA and completed inv()
           across all algebras; built out generic rigid-body dynamics
           (static_/kinematic_/dynamic_system{2,3}dp) with additive force elements
           (joint spring/damper, time-varying applied wrench, grounded spatial spring)
           and an opt-in closed_loop_system KKT constraint layer, all bound into
           ga_py and ga_lua; reorganized the project down to a self-contained generic
           library -- the Qt6 viewer and the application-specific test/model
           bundles were moved out of this repository, the generic dynamic_system{2,3}dp
           gained a protected extra_wrenches() subclass seam so an application can inject
           its own force elements without the base knowing about them, and the ga_docu
           documentation now ships as a generated PDF (its LaTeX sources are maintained
           outside this repository)
- 2026/09: BEHAVIOUR CHANGE -- operator== on every value type is now a RELATIVE
           comparison: equal when each component differs by less than
           eps * max(|a|_inf, |b|_inf, 1), i.e. a budget of 5-10 ulps of the operand at
           any magnitude, with the old absolute window retained below unit scale. The
           previous absolute tolerance was right near magnitude 1 and stricter than the
           floating-point grid above ~1e3 -- at an earth radius one ulp (9.3e-10) exceeded
           the fixed 1.11e-15 window by six orders of magnitude, so values agreeing to
           every digit a double carries compared unequal and is_close had to be reached
           for. The rule now lives in ONE place (detail::coeffs_equal), which the twelve
           operator== bodies call; is_close keeps the larger eps_congruent budget for
           independently computed operands. Note operator== is deliberately approximate
           and therefore NOT transitive: it must not be used to sort, to std::unique or
           to key an associative container. Also in this release: projections,
           reflections, antiprojections and inversions in PGA divide out the scale of the
           object projected/reflected on (they were quadratic in it and silently correct
           only for a unitized target); try_unitize() added for objects that may be
           ideal; central_antiproj2dp/3dp completes the projection family; and sta4ds
           nrm() returns a plain scalar like ega's, leaving the strong norm types to PGA
           where the metric separates bulk from weight
- 2026/09: the projection/rejection grade matrix completed, and its gates written:
           new overloads reject_from(BiVec3dp, TriVec3dp) (pga3dp), project_onto and
           reject_from for Vec4ds and BiVec4ds onto TriVec4ds (sta4ds), and
           project_onto/reject_from(BiVec3d, BiVec3d) (ega3d). Every algebra now has a
           projection-contract test case pinning both scaling contracts, containment,
           idempotence and the degenerate target. Two lessons are recorded at the
           functions: the sign of the textbook (A contracted into B) * inv(B) depends on
           BOTH grades here, because dot folds in the source's reversion -- an even-grade
           source needs the target as B / nrm_sq(B), not inv(B), and only idempotence sees
           the difference; and rejection = source - projection stays a BLADE only when the
           target is a hyperplane, which is why no bivector-onto-bivector pair exists in a
           4d algebra (the remainder there is non-simple, so it is no 2-plane at all)
- 2026/09: projections are computed from ONE named operation per algebra, and carry no
           grade-dependent sign anywhere. Where the grades differ it is the orthogonal
           projection rwdg(b, expansion(a,b)) over the target's squared norm --
           ortho_proj2dp/3dp (PGA, weight dual + weight_nrm_sq) and the new
           ortho_proj3d (EGA3D) / ortho_proj4ds (STA), with the metric dual + nrm_sq;
           where the grades are equal it is dot(a,b)/nrm_sq(b)*b, the only form a
           projection onto a blade of its own grade can take. project_onto delegates
           instead of repeating the expression, in PGA too. The classical
           (b >> a) * inv(b) is gone: its contraction and its inverse each fold in a
           reversion, so an even-grade source needed a (-1)^(j(j-1)/2) that nothing but
           idempotence could catch. Where the grades differ the new form is also faster
           (sta4ds vec->bivec 8.40 -> 2.07 ns, vec->trivec 5.03 -> 1.67 ns), because the
           classical product materialised grades that were then discarded. The ga_lua
           prelude, which still emulated the pre-try_unitize semantics, now matches the
           library
- 2026/09: is_simple() is RELATIVE, and guards the projections. It compared |B ^ B|
           against an absolute 5*machine_eps, a quantity that grows with the components
           squared, so it answered FALSE for genuine 2-blades from unit scale upward (its
           tests used basis blades, where the cancellation is bit-exact, and never saw
           it). It now divides by the sum of the squared components, which is exactly
           scale-invariant and stays defined where nrm_sq cannot serve as a gauge -- a
           Lorentzian metric has null blades, a degenerate one ideal lines.
           is_simple(BiVec3dp) is new (the Pluecker condition: is this bivector a line at
           all?). Both take the tolerance as a parameter defaulting to eps_congruent.
           Projections and reflections that can be given a bivector target in a 4d algebra
           now CHECK it and throw on a non-blade under _HD_GA_EXTENDED_TEST_BLADE_TARGET
           (on by default, droppable per target via ga_blade_guard(<target> OFF)); the
           check selects itself with if constexpr (requires { is_simple(b); }), so exactly
           the types that can fail are tested. It immediately found two places in this
           repository's own tests that projected onto a non-blade, one of them l1 + l2 --
           the sum of two lines, which in 4d is not a line
- 2026/09: the STA causal predicates are RELATIVE too. is_lightlike() tested B^2 == 0
           EXACTLY, so a null vector that was computed rather than typed was never null:
           boosting one leaves B^2 ~ -9e-16 at unit scale and ~ -4e-6 at 1e6, and since
           is_timelike/is_spacelike used strict comparisons, such a vector was classified
           by the SIGN OF THE ROUNDING. normalize() then divided it by nrm ~ 1e-8 instead
           of returning it unchanged, amplifying that rounding by eight orders of
           magnitude. The three now share one threshold relative to the sum of the squared
           components (the same gauge is_simple uses, and exactly scale-invariant since
           both quantities are quadratic), stay mutually exclusive and exhaustive, and
           take the tolerance as a parameter defaulting to eps_congruent. log() and
           angle() / rapidity() improve with them: a near-null plane now takes the branch
           whose limit is exp(B) = 1 + B, and near-null arguments are rejected instead of
           returning an ill-conditioned number
- 2026/09: inv() and normalize() judge degeneracy RELATIVELY (EGA, PGA, STA; CGA still to
           follow). The guard compared the divisor against an absolute 5*machine_eps, so
           inv() of a vector of magnitude 1e-9 threw although its inverse (7e7) is
           perfectly representable -- the same geometry in millimetres instead of metres
           failed, while a genuinely degenerate blade at a large scale passed. What makes
           a blade non-invertible is degeneracy, not size, so the divisor is now judged
           against the sum of the squared coefficients of the object it came from, raised
           to the divisor's own degree (the multivector inverses divide by a
           Hitzer-Sangwine determinant, which is quartic). In PGA that sum is
           bulk_nrm_sq + weight_nrm_sq, in STA the new coeff_sq(); in a EUCLIDEAN metric
           the divisor IS that sum, so the question collapses to "is the object zero" and
           the sites there call check_nonzero() instead -- the honest spelling. Two
           helpers replace check_normalization at 42 sites; a gate pins inv()/normalize()
           down to magnitude 1e-12 and the zero blade to a throw
- 2026/09: the same relative judgement in CGA, which completes it. The 17 inv() sites take
           the metric square against coeff_sq() raised to the divisor's degree, and with
           the null-pair metric that is not a formality: nrm_sq vanishes on every embedded
           point by construction, so the gauge has to be metric-free to exist at all.
           Three more families followed, each a RATIO of two quadratic forms and therefore
           independent of the object's weight -- radius_sq() (the antidot square over the
           round weight), is_flat()/is_round() (whose total_sq IS the sum of the squared
           coefficients, since the four-part split partitions them, so it is the gauge and
           not a threshold) and the generator lines of get_rotation() and
           get_loxodromic(). A sphere whose coefficients were of order 1e-8 has a round
           weight of order 1e-16 and used to throw out of all of them. The weight
           divisions -- unitize() and the dehomogenizing position() -- KEEP an absolute
           floor, because no relative form separates "ideal" from "very far away": the two
           differ exactly by that absolute comparison of weight against bulk. They now say
           so by calling the same check_unitization() PGA uses for them, which leaves
           check_normalization() to the two divisors that are genuinely dimensionless (a
           series value and a dilation factor, both of natural scale one). Gated per
           algebra by "a small object is not a degenerate one", at weights down to 1e-12;
           each conversion was falsified by restoring its absolute form
- 2026/09: ga_py links the ga library instead of reproducing its include path. The two
           development safety gates are INTERFACE compile definitions on the ga target, so
           every consumer gets them by linking it -- but ga_py added the include directory
           by hand and hard-coded _HD_GA_EXTENDED_TEST_DIV_BY_ZERO, so the extension had
           the division guard (copied) and not the blade-target guard (added later,
           never copied): project_onto(point, non-blade) returned a plausible wrong answer
           through the bindings while throwing in C++. ga_lua and ga_view link the target
           and always had both. Copying a definition rather than linking the library that
           owns it opts a target out of every gate added afterwards
- 2026/09: the generic templates the binding generator cannot type are bound from C++ now,
           not reimplemented in Python. The projection / antiprojection families (and
           ega.ortho_proj3d, sta.ortho_proj4ds), dist2dp / dist3dp, try_unitize and
           is_simple's defaulted tolerance live in the new hand-written
           ga_py/src/bindings_projections.cpp, one def per meaningful grade pair. The
           scanner DOES see these functions -- FUNCTION_TEMPLATE is collected beside
           FUNCTION_DECL -- but they are declared as decltype(auto) f(arg1&&, arg2&&), so
           the manifest carries no types to map onto a nanobind signature; project_onto is
           a template too and binds automatically, because it names concrete types. The
           Python reimplementation they replace had drifted twice, both measured through
           the bindings: it divided by an IDEAL target's weight (the library divides only
           above safe_epsilon^2, so the two stood 1e32 apart at weight_nrm_sq = 1e-32),
           and it had no blade-target guard, returning a plausible wrong answer where the
           library throws. try_unitize was not bound at all -- its bool* out-parameter has
           no nanobind mapping -- although the projection documentation tells the reader
           to call it; in Python it returns (object, unitized), the derived primitives
           registered first so a point3dp comes back a point3dp. central_antiproj2dp /
           central_antiproj3dp were never exposed, leaving the four-way family three
           quarters bound. The bulk/weight contractions and expansions deliberately STAY
           Python forwarders: each is a one-line composition of bound primitives,
           identical to its C++ body, and covers every grade pair the composed rwdg / wdg
           / dual accept -- enumerating pairs would narrow them. Five new gates, each
           falsified by restoring the behaviour it pins
- 2026/09: the Euclidean content of a pga object, as nine functions. pga lifts an ega
           quantity by a constructor -- Vec3dp(v, 1) for a point, Vec3dp(v, 0) for a
           direction, BiVec3dp(dir, mom) for a line -- but nothing took the return trip:
           no function in the pga surface returned an ega type, so a consumer computing in
           ega, binding the result to a location in pga and needing it back in ega again
           wrote the components out by hand. to_vec{2,3}d, to_bivec3d, to_scalar3d and
           to_pscalar2d close it, and WHICH part comes back is decided by the RETURN TYPE,
           since every grade has at most one part of each Euclidean type: a line's
           to_vec3d is its direction and its to_bivec3d its moment, a plane's to_vec3d its
           normal and its to_scalar3d its offset (dot(normal, p) + offset == 0 on the
           plane, so z == 3 is TriVec3dp(0, 0, 1, -3)). Nothing is normalized.
           to_vec3d(Vec3dp) is the one that has to ask about the weight, and it asks the
           way try_unitize does: an object carrying weight is a POINT, so its position is
           the bulk over the weight, and an ideal one is a DIRECTION, which is the bulk.
           The test is detail::has_weight on the squared weight, shared with try_unitize
           and detail::by_weight_sq, so the three never disagree -- and the result matches
           unitize() bit for bit. The hand-written vec3d(P.x / P.w, ...) it replaces
           differs by up to 4 ulp (three divisions against a reciprocal and three
           multiplications) and yields Vec3d(inf, inf, inf) for an ideal point. Bound in
           both wrappers: ga_py by a plain regeneration (the signatures are concrete, so
           the generator takes them), ga_lua by hand beside att() and sup()
- 2026/09: BiVec3dp's two-vector ctor is TAG-CONSTRAINED. It takes (Vec3d, BiVec3d) -- a
           direction and a moment -- where it used to accept any two Vec3_t through a
           same-tag and a different-tag overload, so line3d(A, B) with two POINTS compiled
           and silently meant "direction = A, moment = B". Measured on
           line3d(point3d{0,0,1}, point3d{1,0,1}): it yields Line3d(0,0,1,1,0,1), for
           which is_simple is FALSE -- not a line at all, since a line's direction and
           moment must be perpendicular, and a projection onto it now throws under the
           blade guard. The join of two points is join(A, B). Exactly one site in the tree
           used the unconstrained form, the library's own test of that ctor
- 2026/09: ga_lua sheds the same projection reimplementation ga_py did. The Lua prelude
           defined ortho_proj*/central_proj*/ortho_antiproj* and dist{2,3}dp itself, with
           the identical drift: `if n ~= 0.0 and n ~= 1.0` divided by an IDEAL target's
           weight where the library divides only above safe_epsilon^2, and there was no
           blade-target check, so projecting onto a non-blade answered instead of
           throwing.
           Both are gone -- the families are bound from C++ by one lambda per meaningful
           grade pair -- and verified against the C++ values: the ideal target now yields
           Vec3dp(1e-32, 2e-32, -3e-16, 1e-32) in Lua as in C++, and the non-blade target
           raises the library's own message. central_antiproj2dp / central_antiproj3dp had
           no binding at all (three quarters of the four-way family), and try_unitize had
           none either -- in Lua it returns (object, unitized), the language's own idiom
           for a result plus a status. Moving them found a defect the untyped prelude had
           hidden: test_forwarders.lua called BOTH antiprojections with their arguments
           the wrong way round (gr(a) < gr(b), which is a projection's precondition), and
           the prelude computed a meaningless expression that the script accepted because
           it only checked that a userdata came back. The typed binding rejects it. What
           stays in the prelude is the other kind: the 16 bulk/weight contractions and
           expansions plus the two sta expansions, each a one-line composition of bound
           primitives that cannot drift and covers every grade pair it accepts

- 2026/09: qp_ls_solve no longer cycles at a degenerate vertex. Two defects, found by a
           consumer counting the levels that ran to max_iter and were returned silently
           (every one a true cycle: 50x the cap, the point unmoved). The step's least
           squares judged A N's rank at A N's own scale, so a direction a level barely
           sees read as full rank and the step was the residual over ~1e-16 (2e16),
           blocked at zero by the row just released, which re-entered -- 660 iterations;
           it is judged at A's scale now, as solve_hierarchy and the working-set
           factorization already did. And a row released on a negative multiplier
           re-blocked the very next step at zero length forever: a lexicographic level's
           objective is singular whenever it has fewer independent rows than the freedom
           left, and the multiplier's sign then means nothing -- such a row is PINNED
           until x moves, which clears every pin (each iteration moves x, grows the
           working set, pins a row, or stops: all finite; a stop with pinned rows is
           optimal on its face). Three measured remedies were rejected first -- a
           dependent active row dropped before the multipliers decide, a relative
           approach test in the ratio step, a proximal damping of the step -- each fixed
           one instance and cycled on another. Gated by the three instances in
           ga_stencil_test (7 / 9 / 7 iterations, feasible, the objective not above its
           start); the 200-system sweep and every consumer unmoved. And exhausting
           max_iter THROWS now, in qp_ls_solve and bvls_solve alike (iters_out is set
           first): a point the iteration did not settle on is not an answer, and
           returned silently it read as one. The gates gained a second, independent
           oracle beside the brute force over every active set: the KKT certificate
           (multipliers on the whole tight set as a non-negative least squares, zero
           residual = optimal), on the 3-variable sweep, on a new sweep of 200 vertices
           built degenerate on purpose (all certified), and on the three instances --
           where one is NOT: a stop with pinned rows 23 % off its optimum, because the
           descent needs two dependent rows released together. Gated as the known
           failing stop it is; the certificate as an in-loop escape was tried and
           cycled elsewhere
- 2026/09: the robotics chapter of ga_docu states the actuation layer that the library
           carries: a joint's range and drive as its own specification with a neutral
           value; the prioritized stack under the drives' box -- the clipped torque solve,
           the acceleration-level box with no joint ever passive, and the least-squares
           program both are special cases of -- with what a saturated joint means in each;
           a demand refused with the number, the reachable part of a task row as the
           torque box's support function; and the solver's guarantees in one paragraph.
           The generic PDF is refreshed.

