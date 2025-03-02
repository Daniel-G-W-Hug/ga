# Hint for application

- **Single include file for the user** (depending on the GA needed):
  - ga/ga_ega.hpp **or**
  - ga/ga_pga.hpp

- In case the library should be **switched between float and double as value type**, this
  can be adjusted in ga/ga_value_t.hpp at ONE single location

- Just for easy user reference following includes are provided at the top level:
  - ga/ga_usr_consts.hpp
  - ga/ga_usr_types.hpp
  - ga/ga_usr_utilities.hpp

    **These three include files are already included automatically by ga/ga_ega.hpp or
    ga/ga_pga.hpp respectively, and thus must not be included by the user again.**

    All other include files that are not intended for direct use by the user of the
    ga library.
