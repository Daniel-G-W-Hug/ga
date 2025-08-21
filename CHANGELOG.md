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
           Reworked pga2dp and implemented pga3dp based on expressions generated using ga_prdxpr
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
