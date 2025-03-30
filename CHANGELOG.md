# Changes

- 2024/02: add first routines for geometric algebra (GA)
- 2024/04: complete geometric products for 2d and 3d EGA incl. ga_test
- 2024/05: add ga_lua for scripting capabilities
- 2024/07: add ga_view as qt-based viewer for 2d applications
- 2024/08: moved EGA (Euclidean GA) in own namespace and created PGA (Projective GA)
- 2024/12: reworked structure of ega2d and ega3d:
                - common base types established
                - missing products added: dot, wdg, rwdg, inner, fat_dot,
                  left- and right complements
                - duality and complements added
                - unified nrm functions for all types
- 2025/01: added ga_prdxpr to automatically generate product expressions for all products
- 2025/02: started to implement parser for string manipulation with ga_prdxpr to
           manipulate and simplify expressions resulting from sandwich products in order
           to optimize performance (still needs completion for automatic simplification)
           Reworked pga2dp and implemented pga3dp based on expressions generated using ga_prdxpr
- 2025/03: implemented pga2dp and pga3dp including tests; reworked everything into  
           definitions as proposed by Lengyel consistently for ega2d, ega3d, pga2dp, pga3dp; startet to implement latex documentation
